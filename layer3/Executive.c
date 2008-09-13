/* 
   A* -------------------------------------------------------------------
   B* This fil econtains source code for the PyMOL computer program
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

#include"Version.h"
#include"main.h"
#include"Base.h"
#include"OOMac.h"
#include"Executive.h"
#include"ObjectMesh.h"
#include"ObjectDist.h"
#include"ObjectSurface.h"
#include"ObjectSlice.h"
#include"ObjectAlignment.h"
#include"ObjectGroup.h"
#include"ListMacros.h"
#include"Ortho.h"
#include"Scene.h"
#include"Selector.h"
#include"Vector.h"
#include"Color.h"
#include"Setting.h"
#include"Matrix.h"
#include"P.h"
#include"PConv.h"
#include"Match.h"
#include"ObjectCGO.h"
#include"Util.h"
#include"Wizard.h"
#include"ScrollBar.h"
#include"Movie.h"
#include"ObjectGadgetRamp.h"
#include"SculptCache.h"
#include"Control.h"
#include"Menu.h"
#include"Map.h"
#include"Editor.h"
#include"RepDot.h"
#include"Seq.h"
#include"Text.h"
#include"PyMOL.h"
#include"PyMOLOptions.h"
#include"Tracker.h"
#include"Word.h"
#include"main.h"
#include"Parse.h"
#include"PlugIOManager.h"

#include"OVContext.h"
#include"OVLexicon.h"
#include"OVOneToOne.h"
#include"OVOneToAny.h"

#define cExecObject 0
#define cExecSelection 1
#define cExecAll 2

#define cTempRectSele "_rect"
#define cLeftButSele "lb"
#define cIndicateSele "indicate"


typedef struct SpecRec { 
  /* NOTE: must zero-init with CALLOC */
  int type;
  WordType  name; /*only used for selections*/
  CObject *obj;  
  struct SpecRec *next;
  int repOn[cRepCnt];
  int visible;

  ObjectNameType group_name;

  /* not pickled */
  int sele_color;
  int hilight; /* 0 = none, 1 = name, 2 = group control (if any) */
  int previous;
  int cand_id;
  struct SpecRec *group;
  int group_member_list_id; 
  int in_scene,is_hidden;
  int in_panel;
  int grid_slot;
} SpecRec; /* specification record (a line in the executive window) */

typedef struct PanelRec {
  SpecRec *spec;
  int nest_level;
  int is_group;
  int is_open;
  struct PanelRec *next;
} PanelRec;

typedef struct {
  int list_id;
  int next;
} ListMember;

struct _CExecutive {
  Block *Block;
  SpecRec *Spec;
  CTracker *Tracker;
  int Width,Height,HowFarDown;
  int ScrollBarActive;
  int NSkip;
  struct CScrollBar *ScrollBar;
  CObject *LastEdited;
  int DragMode;
  int Pressed,Over,LastOver,OldVisibility,ToggleMode,PressedWhat,OverWhat;
  SpecRec *LastChanged,*LastZoomed,*RecoverPressed;
  int ReorderFlag;
  OrthoLineType ReorderLog;
  int oldPX,oldPY,oldWidth,oldHeight,sizeFlag;
  int all_names_list_id, all_obj_list_id, all_sel_list_id;
  OVLexicon *Lex;
  OVOneToOne *Key;
  int ValidGroups;
  int ValidSceneMembers;
  int ValidGridSlots;
  PanelRec *Panel;
  int ValidPanel;
  int CaptureFlag;
};

/* routines that still need to be updated for Tracker list iteration

Low priority:

ExecutiveSculptIterate
ExecutiveSculptActivate
ExecutiveSculptDeactivate
ExecutiveRenameObjectAtoms
ExecutiveSpheroid

*/

static void ExecutiveSpecEnable(PyMOLGlobals *G, SpecRec *rec, int parents, int log);
static void ExecutiveToggleAllRepVisib(PyMOLGlobals *G,int rep);
static void ExecutiveSetAllRepVisib(PyMOLGlobals *G,int rep,int state);
static SpecRec *ExecutiveFindSpec(PyMOLGlobals *G,char *name); 
static int ExecutiveDrag(Block *block,int x,int y,int mod);
static void ExecutiveSpecSetVisibility(PyMOLGlobals *G,SpecRec *rec,
                                       int new_vis,int mod,int parents);
static int ExecutiveSetObjectMatrix2(PyMOLGlobals *G,CObject *obj,int state,double *matrix);
static int ExecutiveGetObjectMatrix2(PyMOLGlobals *G,CObject *obj,int state,double **matrix, int incl_ttt);
int ExecutiveTransformObjectSelection2(PyMOLGlobals *G,CObject *obj,int state,
                                       char *s1,int log,float *matrix,
                                       int homogenous,int global);

int ExecutiveIsosurfaceEtc(PyMOLGlobals *G, 
                           char *surf_name, char *map_name, float lvl, 
                           char *sele, float fbuf, int state, 
                           float carve, int map_state, int side,
                           int quiet, int surf_mode, int box_mode) 
{
  int c;
  OrthoLineType s1;
  CObject *obj=NULL,*mObj,*origObj;
  ObjectMap *mapObj;
  float mn[3] = { 0,0,0};
  float mx[3] = { 15,15,15};
  float *vert_vla = NULL;
  int ok = false;
  ObjectMapState *ms;
  int multi=false;
  /* box_mode 0 = all, 1 = sele + buffer, 2 = vector, 3 = testing */

  
  origObj=ExecutiveFindObjectByName(G,surf_name);  
  if(origObj) {
    if(origObj->type!=cObjectSurface) {
      ExecutiveDelete(G,surf_name);
      origObj=NULL;
    }
  }
  
  mObj=ExecutiveFindObjectByName(G,map_name);  
  if(mObj) {
    if(mObj->type!=cObjectMap)
      mObj=NULL;
  }
  if(mObj) {
    mapObj = (ObjectMap*)mObj;
    if(state==-1) {
      multi=true;
      state=0;
      map_state=0;
    } else if(state==-2) { /* current state */
      state=SceneGetState(G);
      if(map_state<0) 
        map_state=state;
    } else if(state==-3) { /* append mode */
      state=0;
      if(origObj)
        if(origObj->fGetNFrame)
          state=origObj->fGetNFrame(origObj);
    } else {
      if(map_state==-1) {
        map_state=0;
        multi=true;
      } else {
        multi=false;
      }
    }
    while(1) {
      if(map_state==-2)
        map_state=SceneGetState(G);
      if(map_state==-3)
        map_state=ObjectMapGetNStates(mapObj)-1;
      ms = ObjectMapStateGetActive(mapObj,map_state);
      if(ms) {
        switch(box_mode) { 
        case 0: /* using map to get extents */
          for(c=0;c<3;c++) {
            mn[c] = ms->Corner[c];
            mx[c] = ms->Corner[3*7+c];
          }
          if(ms->State.Matrix) {
            transform44d3f(ms->State.Matrix,mn,mn);
            transform44d3f(ms->State.Matrix,mx,mx);
            {
              float tmp;
              int a;
              for(a=0;a<3;a++)
                if(mn[a]>mx[a]) { tmp=mn[a];mn[a]=mx[a];mx[a]=tmp; }
            }
          }
          carve = 0.0F;
          break;
        case 1: /* using selection to get extents */
          ok = (SelectorGetTmp(G,sele,s1)>=0);
          ExecutiveGetExtent(G,s1,mn,mx,false,-1,false); /* TODO state */
          if(carve!=0.0F) {
            vert_vla = ExecutiveGetVertexVLA(G,s1,state);
            if(fbuf<=R_SMALL4)
              fbuf = fabs(carve);
          }
          SelectorFreeTmp(G,s1);
          for(c=0;c<3;c++) {
            mn[c]-=fbuf;
            mx[c]+=fbuf;
          }
          break;
        }
        PRINTFB(G,FB_CCmd,FB_Blather)
          " Isosurface: buffer %8.3f carve %8.3f\n",fbuf,carve
          ENDFB(G);
        obj=(CObject*)ObjectSurfaceFromBox(G,(ObjectSurface*)origObj,mapObj,map_state,
                                           state,mn,mx,lvl,surf_mode,
                                           carve,vert_vla,side,quiet);
        /* copy the map's TTT */
        ExecutiveMatrixCopy2(G, 
                             mObj, obj, 1, 1, 
                             -1,-1, false, 0, quiet);

        if(!origObj) {
          ObjectSetName(obj,surf_name);
          ExecutiveManageObject(G,(CObject*)obj,-1,quiet);
        }
        if(SettingGet(G,cSetting_isomesh_auto_state))
          if(obj) ObjectGotoState((ObjectMolecule*)obj,state);
        if(!quiet) {
          PRINTFB(G,FB_ObjectSurface,FB_Actions)
            " Isosurface: created \"%s\", setting level to %5.3f\n",surf_name,lvl
            ENDFB(G);
        }
      } else if(!multi) {
        PRINTFB(G,FB_ObjectMesh,FB_Warnings)
          " Isosurface-Warning: state %d not present in map \"%s\".\n",map_state+1,map_name
          ENDFB(G);
        ok = false;
      }
      if(multi) {
        origObj = obj;
        map_state++;
        state++;
        if(map_state>=mapObj->NState)
          break;
      } else {
        break;
      }
    }
  } else {
    PRINTFB(G,FB_ObjectSurface,FB_Errors)
      " Isosurface: Map or brick object \"%s\" not found.\n",map_name
      ENDFB(G);
    ok = false;
  }
  return ok;
}

int ExecutiveIsomeshEtc(PyMOLGlobals *G, 
                        char *mesh_name, char *map_name, float lvl, 
                        char *sele, float fbuf, int state, 
                        float carve, int map_state, int quiet,
                        int mesh_mode, int box_mode, float alt_lvl) 
{
  int ok=true;
  CObject *obj=NULL,*mObj,*origObj;
  ObjectMap *mapObj;
  float mn[3] = { 0,0,0};
  float mx[3] = { 15,15,15};
  float *vert_vla = NULL;
  int multi=false;
  ObjectMapState *ms;
  OrthoLineType s1;
  ObjectMolecule *sele_obj = NULL;

  origObj=ExecutiveFindObjectByName(G,mesh_name);  
  if(origObj) {
    if(origObj->type!=cObjectMesh) {
      ExecutiveDelete(G,mesh_name);
      origObj=NULL;
    }
  }
  
  mObj=ExecutiveFindObjectByName(G,map_name);  
  if(mObj) {
    if(mObj->type!=cObjectMap)
      mObj=NULL;
  }
  if(mObj) {
    mapObj = (ObjectMap*)mObj;
    if(state==-1) {
      multi=true;
      state=0;
      map_state=0;
    } else if(state==-2) {
      state=SceneGetState(G);
      if(map_state<0) 
        map_state=state;
    } else if(state==-3) { /* append mode */
      state=0;
      if(origObj)
        if(origObj->fGetNFrame)
          state=origObj->fGetNFrame(origObj);
    } else {
      if(map_state==-1) {
        map_state=0;
        multi=true;
      } else {
        multi=false;
      }
    }
    while(1) {
      if(map_state==-2)
        map_state=SceneGetState(G);
      if(map_state==-3)
        map_state=ObjectMapGetNStates(mapObj)-1;
      ms = ObjectMapStateGetActive(mapObj,map_state);
      if(ms) {
        switch(box_mode) {
        case 0: /* do the whole map */
          {
            int c;
            for(c=0;c<3;c++) {
              mn[c] = ms->Corner[c];
              mx[c] = ms->Corner[3*7+c];
            }
          }
          if(ms->State.Matrix) {
            transform44d3f(ms->State.Matrix,mn,mn);
            transform44d3f(ms->State.Matrix,mx,mx);
            {
              float tmp;
              int a;
              for(a=0;a<3;a++)
                if(mn[a]>mx[a]) { tmp=mn[a];mn[a]=mx[a];mx[a]=tmp; }
            }
          }
          carve = -0.0; /* impossible */
          break;
        case 1: /* just do area around selection */
          ok = (SelectorGetTmp(G,sele,s1)>=0);
          if(ok) {
            int sele1 = SelectorIndexByName(G,s1);
            if(sele1>=0) sele_obj = SelectorGetSingleObjectMolecule(G,sele1);
          }
            
          ExecutiveGetExtent(G,s1,mn,mx,false,-1,false); /* TODO state */
          if(carve!=0.0) {
            vert_vla = ExecutiveGetVertexVLA(G,s1,state);
            if(fbuf<=R_SMALL4)
              fbuf = fabs(carve);
          }
          SelectorFreeTmp(G,s1);
          {
            int c;
            for(c=0;c<3;c++) {
              mn[c]-=fbuf;
              mx[c]+=fbuf;
            }
          }
          break;
        }
        PRINTFB(G,FB_CCmd,FB_Blather)
          " Isomesh: buffer %8.3f carve %8.3f \n",fbuf,carve
          ENDFB(G);
        if(sele_obj && SettingGet_b(G,NULL,sele_obj->Obj.Setting,cSetting_map_auto_expand_sym) &&
           (sele_obj->Symmetry) && ObjectMapValidXtal(mapObj,state)) { 
          obj=(CObject*)ObjectMeshFromXtalSym(G,(ObjectMesh*)origObj,mapObj,
                                           sele_obj->Symmetry,
                                           map_state,state,mn,mx,lvl,mesh_mode,
                                           carve,vert_vla,alt_lvl,quiet);
        } else {
          obj = NULL;
        }
        if(!obj) {
          obj=(CObject*)ObjectMeshFromBox(G,(ObjectMesh*)origObj,mapObj,
                                          map_state,state,mn,mx,lvl,mesh_mode,
                                          carve,vert_vla,alt_lvl,quiet);
        }
        /* copy the map's TTT */
        ExecutiveMatrixCopy2(G, 
                             mObj, obj, 1, 1, 
                             -1,-1, false, 0, quiet);
        
        if(!origObj) {
          ObjectSetName(obj,mesh_name);
          ExecutiveManageObject(G,(CObject*)obj,false,quiet);
        }          
        
        if(SettingGet(G,cSetting_isomesh_auto_state))
          if(obj) ObjectGotoState((ObjectMolecule*)obj,state);
        if(!quiet) {
          if(mesh_mode!=3) {
            PRINTFB(G,FB_ObjectMesh,FB_Actions)
              " Isomesh: created \"%s\", setting level to %5.3f\n",mesh_name,lvl
              ENDFB(G);
          } else {
            PRINTFB(G,FB_ObjectMesh,FB_Actions)
              " Gradient: created \"%s\"\n",mesh_name
              ENDFB(G);
          }
        }
      } else if(!multi) {
        PRINTFB(G,FB_ObjectMesh,FB_Warnings)
          " Isomesh-Warning: state %d not present in map \"%s\".\n",map_state+1,map_name
          ENDFB(G);
        ok = false;
      }
      if(multi) {
        origObj = obj;
        map_state++;
        state++;
        if(map_state>=mapObj->NState)
          break;
      } else {
        break;
      }
    }
  } else {
    PRINTFB(G,FB_ObjectMesh,FB_Errors)
      " Isomesh: Map or brick object \"%s\" not found.\n",map_name
      ENDFB(G);
    ok = false;
  }
  return ok;
}

int ExecutivePseudoatom(PyMOLGlobals *G, char *object_name, char *sele,
                        char *name, char *resn, char *resi, char *chain,
                        char *segi, char *elem, float vdw, int hetatm,
                        float b, float q, char *label, float *pos, int color, 
                        int state, int mode,  int quiet)
{
  int ok = true;
  
  ObjectMolecule *obj = ExecutiveFindObjectMoleculeByName(G,object_name);
  int is_new = false;
  int sele_index = -1;
  float local_pos[3];

  if(sele && sele[0]) {
    if(WordMatch(G,cKeywordCenter,sele,1)<0) {
      sele = NULL;
      SceneGetPos(G,local_pos);
      pos = local_pos;
    } else if(WordMatch(G,cKeywordOrigin,sele,1)<0) {
      sele = NULL;
      SceneOriginGet(G,local_pos);
      pos = local_pos;
    }
  }

  if(sele && sele[0]) {
    sele_index=SelectorIndexByName(G,sele);    
    if(sele_index<0) {
      ok = false;
      PRINTFB(G,FB_Executive,FB_Errors)
        " Pseudoatom-Error: invalid selection\n"
        ENDFB(G);
    }
  }
  if(ok) {
    if(!obj) {
      /* new object */
      is_new = true;
      obj = ObjectMoleculeNew(G,false);
      ObjectSetName(&obj->Obj,object_name);
      if(!obj)
        ok=false;
    }
  }

  if(ok) {
    if(ObjectMoleculeAddPseudoatom(obj,sele_index, name, resn, resi, chain,
                                   segi, elem, vdw, hetatm, b, q, label, pos, color,
                                   state, mode, quiet)) {
      if(is_new) {
        ExecutiveDelete(G,object_name); /* just in case */
        ExecutiveManageObject(G,&obj->Obj,false,true);
      } else {
        ExecutiveUpdateObjectSelection(G,&obj->Obj);
      }
    }
  }
  return ok;
}


static void ExecutiveInvalidateGridSlots(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  I->ValidGridSlots = false;
}

static void ExecutiveUpdateGridSlots(PyMOLGlobals *G, int force)
{
  register CExecutive *I = G->Executive;
  int grid_slot_count = 0;
  int grid_by_group = 1; /* grid slots are inherited this many levels */

  ExecutiveUpdateGroups(G,false);
  if(force || (!I->ValidGridSlots)) {
    CTracker *I_Tracker= I->Tracker;
    I->ValidGridSlots = true;
    {
      SpecRec *rec=NULL;
      while(ListIterate(I->Spec,rec,next)) {
        rec->grid_slot = 0;
        if(rec->type==cExecObject) {
          /* make sure every object (potentially) needing a grid slot gets one */
          switch(rec->obj->type) {
          case cObjectMolecule:
          case cObjectMap:
          case cObjectMesh:
          case cObjectMeasurement:
          case cObjectCallback:
          case cObjectCGO:
          case cObjectSurface:
          case cObjectSlice:
          case cObjectGadget:
          case cObjectGroup:
            if(!rec->grid_slot)
              rec->grid_slot = ++grid_slot_count; 
            break;
          }
        }
      }
    }
    
    if(grid_by_group) {
      SpecRec *rec=NULL,*group_rec = NULL;
      while(ListIterate(I->Spec,rec,next)) {
        OVreturn_word result;
        if( OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,rec->group_name)))) {
          if( OVreturn_IS_OK( (result = OVOneToOne_GetForward(I->Key, result.word)))) { 
            if(TrackerGetCandRef(I_Tracker, result.word, (TrackerRef**)&group_rec)) {
              register int grid_slot_group_depth = grid_by_group;
              { 
                SpecRec *check_rec = group_rec;
                while(check_rec && grid_slot_group_depth) {
                  if(grid_slot_group_depth==1)
                    rec->grid_slot = check_rec->grid_slot;
                  if(check_rec == rec) { /* cycle */
                    break;
                  } else {
                    check_rec = check_rec->group;
                    grid_slot_group_depth--;
                  }
                }
              }
            }
          }
        }
      }
    }
      
    {
      SpecRec *rec=NULL;
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject) {
          int obj_slot = SettingGet_i(G,rec->obj->Setting,NULL,cSetting_grid_slot);
          if(obj_slot == -1) {
            rec->obj->grid_slot = rec->grid_slot;
          } else
            rec->obj->grid_slot = obj_slot;
        }
      }
    }
  }
}

static void ExecutiveInvalidatePanelList(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  if(I->ValidPanel) {
    if(I->Panel) 
      ListFree(I->Panel,next,PanelRec);
    I->ValidPanel = false;
    ExecutiveInvalidateGridSlots(G);
  }
}


static PanelRec *PanelListGroup(PyMOLGlobals *G, PanelRec *panel, SpecRec *group,
                                int level,int hide_underscore)
{
  register CExecutive *I = G->Executive;
  PanelRec *result = NULL;
  SpecRec *rec = NULL;
  /* set up recursion prevention */
  while(ListIterate(I->Spec,rec,next)) {
    rec->in_panel = false;
  }
  while(ListIterate(I->Spec,rec,next)) { /* add all members which belong to this group */

    if((rec->name[0]!='_')||(!hide_underscore)) { /* not hidden */
      if((rec->group == group)&&(!rec->in_panel)) {
        int group_name_len = strlen(rec->group_name);
        if((!hide_underscore)||!
           ((strncmp(rec->name,rec->group_name,group_name_len)==0) && /* named with proper group prefix */
            (rec->name[group_name_len]=='.') &&
            (rec->name[group_name_len+1]=='_'))) { /* and not hidden inside group */
          
          PanelRec *new_panel = NULL;
          ListElemCalloc(G,new_panel,PanelRec);
          if(panel) 
            panel->next = new_panel;
          else
            result = new_panel;
          panel = new_panel;
          panel->spec = rec;
          panel->nest_level = level;
          if(!level) rec->group_name[0] = 0; /* force open any cycles which have been created...*/
          rec->in_panel = true;
          if((rec->type == cExecObject) && 
             (rec->obj->type == cObjectGroup)) {
            ObjectGroup *obj_group = (ObjectGroup*)rec->obj;
            panel->is_group = true;
            if(obj_group->OpenOrClosed) {
              panel->is_open = true;
              panel = PanelListGroup(G,panel,rec,level+1,hide_underscore);
            }
          }
        }
      }
    }
  }
  if(!result)
    result = panel;
  return result;
}

static void ExecutiveUpdatePanelList(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  int hide_underscore = SettingGetGlobal_b(G,cSetting_hide_underscore_names);
  if(!I->ValidPanel) {
    /* brute-force & inefficient -- need to optimize algorithm */
    I->Panel = PanelListGroup(G,NULL,NULL,0,hide_underscore);
    I->ValidPanel = true;
  }
}

static void ExecutiveInvalidateSceneMembers(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  I->ValidSceneMembers=false;
}

static void ExecutiveUpdateSceneMembers(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  ExecutiveUpdateGroups(G,false);
  ExecutiveUpdateGridSlots(G,false);
  if(!I->ValidSceneMembers) {
    SpecRec *rec=NULL;
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        int visible = rec->visible;
        SpecRec *group_rec = rec->group;
        while(visible && group_rec) { /* visibility is a group issue... */
          if(!group_rec->visible)
            visible=false;
          else
            group_rec = group_rec->group;
        }
        if(rec->in_scene && !visible) {
          rec->in_scene = SceneObjectDel(G,rec->obj);
        } else if(visible && !rec->in_scene) {
          rec->in_scene = SceneObjectAdd(G,rec->obj);
        }
      }
    }
    I->ValidSceneMembers = true;
  }
}

void ExecutiveInvalidateGroups(PyMOLGlobals *G,int force)
{
  register CExecutive *I = G->Executive;
  if(force || I->ValidGroups) {
    CTracker *I_Tracker= I->Tracker;
    SpecRec *rec=NULL;
    while(ListIterate(I->Spec,rec,next)) {
      rec->group = NULL;
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectGroup) {
          int list_id = rec->group_member_list_id;
          if(list_id) 
            TrackerDelList(I_Tracker, rec->group_member_list_id);
          rec->group_member_list_id = 0; /* not a list */
        }
    }
    I->ValidGroups = false;
    ExecutiveInvalidateSceneMembers(G); 
    ExecutiveInvalidatePanelList(G);
  }
  /* any changes to group structure means that we need to check scene
     members */
}


void ExecutiveUpdateGroups(PyMOLGlobals *G,int force)
{
  register CExecutive *I = G->Executive;

  if(force || (!I->ValidGroups)) {
    CTracker *I_Tracker= I->Tracker;

    /* first, get rid of existing group lists */
    
    if(force || I->ValidGroups) ExecutiveInvalidateGroups(G,true);

    /* create empty lists for each group (also init grid_slot) */
    
    {
      SpecRec *rec=NULL;
      while(ListIterate(I->Spec,rec,next)) {
        rec->group = NULL;
        if(rec->type==cExecObject) {
          if(rec->obj->type==cObjectGroup) {
            rec->group_member_list_id = TrackerNewList(I_Tracker,NULL);
          }
        }
      }
    }
              
    /* iterate through and populate groups lists with their members */

    {
      SpecRec *rec=NULL,*group_rec = NULL;
      while(ListIterate(I->Spec,rec,next)) {
        OVreturn_word result;
        if( OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,rec->group_name)))) {
          if( OVreturn_IS_OK( (result = OVOneToOne_GetForward(I->Key, result.word)))) { 
            if(TrackerGetCandRef(I_Tracker, result.word, (TrackerRef**)&group_rec)) {
              int cycle = false;
              { /* don't close infinite loops */
                SpecRec *check_rec = group_rec;
                while(check_rec) {
                  if(check_rec == rec) {
                    cycle=true;
                    break;
                  } else {
                    check_rec = check_rec->group;
                  }
                }
              }
              if(!cycle) {
                rec->group = group_rec;
                TrackerLink(I_Tracker,rec->cand_id,group_rec->group_member_list_id,1);
              }
            }
          }
        }
      }
    }

    {
      int hide_underscore = SettingGetGlobal_b(G,cSetting_hide_underscore_names);
      if(hide_underscore) {
        SpecRec *rec=NULL;
        while(ListIterate(I->Spec,rec,next)) {
          rec->is_hidden=false;
          if(rec->name[0]=='_')
            rec->is_hidden=true;
          else if(rec->group) {
            int group_name_len = strlen(rec->group_name);
            if(rec->group->is_hidden)
              rec->is_hidden=true;
            else if((strncmp(rec->name,rec->group_name,group_name_len)==0) && 
                    (rec->name[group_name_len]=='.') &&
                    (rec->name[group_name_len+1]=='_'))
              rec->is_hidden=true;
          }
        }
        { /* sub-optimal propagation of hidden status to group members */
          int repeat_flag=true;
          while(repeat_flag) {
            repeat_flag=false;
            while(ListIterate(I->Spec,rec,next)) {
              if(rec->group&&(!rec->is_hidden)) {
                if(rec->group->is_hidden) {
                  rec->is_hidden = true;
                  repeat_flag = true;
                }
              }
            }
          }
        }
      }
    }

    /* note that it is possible to have infinite loops -- these must be
       allowed for later in the group expansion routine(s) */
    I->ValidGroups = true;
    ExecutiveInvalidatePanelList(G);
  }
}

static int ExecutiveGetObjectParentList(PyMOLGlobals *G, SpecRec *child)
{
  int list_id = 0;
  ExecutiveUpdateGroups(G,false);
  {
    CExecutive *I = G->Executive;
    CTracker *I_Tracker= I->Tracker;
    int priority = 1; /* generations removed from child */
    int repeat_flag = true;
    SpecRec *group_rec = NULL;

    list_id = TrackerNewList(I_Tracker, NULL);
    while(child && child->group && repeat_flag) {
      OVreturn_word result;
      repeat_flag = false;
      if( OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,child->group_name)))) {
        if( OVreturn_IS_OK( (result = OVOneToOne_GetForward(I->Key, result.word)))) { 
          if(TrackerGetCandRef(I_Tracker, result.word, (TrackerRef**)&group_rec)) {
            if(TrackerLink(I_Tracker, result.word,  list_id, priority++)) { /* checking this prevents infinite loops */
              if(group_rec->group) {
                repeat_flag=true;
                child = group_rec;
              }
            }
          }
        }
      }
    }
  }
  return list_id;
}

int ExecutiveVdwFit(PyMOLGlobals *G,char *s1,int state1,char *s2,int state2,float buffer, int quiet)
{
  int sele1=SelectorIndexByName(G,s1);
  int sele2=SelectorIndexByName(G,s2);
  int ok=true;

  if((sele1>=0)&&(sele2>=0)) {
    ok = SelectorVdwFit(G,sele1,state1,sele2,state2,buffer,quiet);
  } else {
    ok =false;
  }
  return ok;
}

static int get_op_cnt(PyMOLGlobals *G)
{
  int result = 5;
  if((SettingGetGlobal_i(G,cSetting_button_mode)==2)&&
     !strcmp(SettingGetGlobal_s(G,cSetting_button_mode_name),"3-Button Motions"))
    result = 6;
  return result;
}

static int ExecutiveAddKey(CExecutive *I, SpecRec *rec)
{
  int ok=false;
  OVreturn_word result;
  if(OVreturn_IS_OK( (result = OVLexicon_GetFromCString(I->Lex,rec->name)))) {
    if(OVreturn_IS_OK(OVOneToOne_Set(I->Key, result.word, rec->cand_id))) {
      ok=true;
    }
  }
  return ok;
}
static int ExecutiveDelKey(CExecutive *I, SpecRec *rec)
{
  int ok=false;
  OVreturn_word result;
  if(OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,rec->name)))) {
    if(OVreturn_IS_OK(OVLexicon_DecRef(I->Lex, result.word)) &&
       OVreturn_IS_OK(OVOneToOne_DelForward(I->Key, result.word))) {
      ok=true;
    }
  }
  return ok;
}

static SpecRec *ExecutiveUnambiguousNameMatch(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  SpecRec *result = NULL;
  SpecRec *rec=NULL;
  int best = 0;
  int wm;
  int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);

  while(ListIterate(I->Spec,rec,next)) {
    wm = WordMatch(G,name,rec->name,ignore_case);
    if(wm<0) { /* exact match, so this is valid */
      result = rec;
      best = wm;
      break;
    } else if ((wm>0)&&(best<wm)) {
      result = rec;
      best = wm;
    } else if ((wm>0)&&(best==wm)) { /* ambiguous match... no good */
      result = NULL;
    }
  }
  return(result);
}

static SpecRec *ExecutiveAnyCaseNameMatch(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  SpecRec *result = NULL;
  SpecRec *rec=NULL;

  int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);
  while(ListIterate(I->Spec,rec,next)) {
    if(WordMatchExact(G,name,rec->name,ignore_case)) {
      result = rec;
      break;
    }
  }
  return(result);
}
void ExecutiveUpdateColorDepends(PyMOLGlobals *G,ObjectMolecule *mol)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      if(rec->obj->type==cObjectGadget) { 
        ObjectGadget *gadget = (ObjectGadget*)rec->obj;
        if(gadget->GadgetType == cGadgetRamp) {
          ObjectGadgetRamp *ramp = (ObjectGadgetRamp*)gadget;
          if(ramp->RampType==cRampMol) {
            if(ramp->Mol == mol) {
              ExecutiveInvalidateRep(G,cKeywordAll,cRepAll,cRepInvColor);
              break;
            }
          }
        }
      }
    }
  }
}

void ExecutiveUpdateCoordDepends(PyMOLGlobals *G,ObjectMolecule *mol)
{ /* nasty, ugly, inefficient hack */

  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      if(rec->obj->type==cObjectGadget) { 
        ObjectGadget *gadget = (ObjectGadget*)rec->obj;
        if(gadget->GadgetType == cGadgetRamp) {
          ObjectGadgetRamp *ramp = (ObjectGadgetRamp*)gadget;
          if(ramp->RampType==cRampMol) {
            if(ramp->Mol == mol) {
              ExecutiveInvalidateRep(G,cKeywordAll,cRepAll,cRepInvColor);
              break;
            }
          }
        }
      }
    }
  }
}

int ExecutiveValidNamePattern(PyMOLGlobals *G,char *name)
{
  int result = false;
  CWordMatcher *matcher;
  CWordMatchOptions options;
  char *wildcard = SettingGetGlobal_s(G,cSetting_wildcard);
      
  WordMatchOptionsConfigNameList(&options, 
                                 *wildcard,
                                 SettingGetGlobal_b(G,cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, false);
  if(matcher) { /* this appears to be a pattern */
    result = true;
    WordMatcherFree(matcher);
  } else if(ExecutiveUnambiguousNameMatch(G,name)) {
    /* this does not appear to be a pattern, so it is an unambiguous partial name? */
    result = true;
  }
  return result;

}

#define cExecNoExpand false
#define cExecExpandGroups true
#define cExecExpandKeepGroups 2

static void ExecutiveExpandGroupsInList(PyMOLGlobals *G,int list_id,int expand_groups)
{
  register CExecutive *I = G->Executive;
  CTracker *I_Tracker= I->Tracker;  
  int new_member_added = true;
  SpecRec *rec;
  ExecutiveUpdateGroups(G,false);
  while(new_member_added) { /* keep adding til we can't add no more */
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int cand_id;
    new_member_added = false;
    if(iter_id) {
      while( (cand_id = TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec)) ) {
        if(rec && (rec->type==cExecObject) && 
           rec->group_member_list_id && (rec->obj->type==cObjectGroup)) {
          int group_iter_id = TrackerNewIter(I_Tracker, 0, rec->group_member_list_id);
          int group_cand_id;  
          SpecRec *group_rec;
          if(group_iter_id) {
            while( (group_cand_id = TrackerIterNextCandInList(I_Tracker, group_iter_id,
                                                              (TrackerRef**)&group_rec)) ) {
              if(group_rec && group_cand_id) {
                if(TrackerLink(I_Tracker, group_cand_id, list_id, 1))
                  new_member_added = true;
              }
            }
            TrackerDelIter(I_Tracker, group_iter_id);
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }
  }
  /* now purge all groups from the expanded list */
  if(expand_groups != cExecExpandKeepGroups) {
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int cand_id;
    while( (cand_id = TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec)) ) {
      if(rec && (rec->type==cExecObject) && (rec->obj->type==cObjectGroup)) {
        TrackerUnlink(I_Tracker, cand_id, list_id);
      }
    }
  }
}


/* DON'T FORGET TO RELEASE LIST WHEN DONE!!! */
static int ExecutiveGetNamesListFromPattern(PyMOLGlobals *G,char *name,
                                            int allow_partial,int expand_groups)
{
  register CExecutive *I = G->Executive;
  int result = 0;
  CWordMatcher *matcher;
  CWordMatchOptions options;
  CTracker *I_Tracker= I->Tracker;
  char *wildcard = SettingGetGlobal_s(G,cSetting_wildcard);
  int iter_id = TrackerNewIter(I_Tracker, 0, I->all_names_list_id);
  int cand_id;
  int group_found = false;
  SpecRec *rec;
      
  WordMatchOptionsConfigNameList(&options, 
                                 *wildcard,
                                 SettingGetGlobal_b(G,cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, false);
  if(matcher) {
    if(iter_id) {
      while( (cand_id = TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec)) ) {
        if(rec && !(rec->type==cExecAll)) {
          if(WordMatcherMatchAlpha(matcher,rec->name)) {
            if((rec->type==cExecObject)&&(rec->obj->type==cObjectGroup))
              group_found=true;
            if(!result)  
              result = TrackerNewList(I_Tracker, NULL);
            if(result) {
              TrackerLink(I_Tracker, cand_id, result, 1);
            }
          }
        }
      }
    }
  } else if( (rec = ExecutiveFindSpec(G,name) ) ) { /* only one name in list */
    if((rec->type==cExecObject) && (rec->obj->type==cObjectGroup))
      group_found=true;
    result = TrackerNewList(I_Tracker, NULL);
    TrackerLink(I_Tracker, rec->cand_id, result, 1);
  } else if( allow_partial && (rec = ExecutiveUnambiguousNameMatch(G,name))) {
    if((rec->type==cExecObject) && (rec->obj->type==cObjectGroup))
      group_found=true;
    result = TrackerNewList(I_Tracker, NULL);
    TrackerLink(I_Tracker, rec->cand_id, result, 1);
  }
  if(matcher) WordMatcherFree(matcher);
  if(iter_id) TrackerDelIter(I->Tracker, iter_id);
  if(group_found && expand_groups) {
    ExecutiveExpandGroupsInList(G,result,expand_groups);
  }
  return result;
}

int ExecutiveGroup(PyMOLGlobals *G,char *name,char *members,int action, int quiet)
{
  int ok=true;
  CExecutive *I = G->Executive;
  CObject *obj = ExecutiveFindObjectByName(G,name);

  if(obj && (obj->type!=cObjectGroup)) {
    if((action!=7)||(members[0])) {
      PRINTFB(G,FB_Executive,FB_Errors)
        " Group-Error: object '%s' is not a group object.", name
        ENDFB(G);
      ok=false;
    }
  } else {
    if((!obj)&&(action==cExecutiveGroupAdd)) {
      obj = (CObject*)ObjectGroupNew(G);
      if(obj) {
        ObjectSetName(obj,name);
        ExecutiveManageObject(G,obj,false,true);
      }
    }
  }
  if((!members[0])&&((action == cExecutiveGroupOpen)||
                     (action == cExecutiveGroupClose)||
                     (action == cExecutiveGroupUngroup)||
                     (action == cExecutiveGroupToggle)||
                     (action == cExecutiveGroupEmpty)||
                     (action == cExecutiveGroupPurge)||
                     (action == cExecutiveGroupExcise))) {
    ExecutiveUpdateGroups(G,false);
    {
      CTracker *I_Tracker= I->Tracker;
      int list_id = ExecutiveGetNamesListFromPattern(G,name,true,false);
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec *rec;
      
      while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
        if(rec) {
          ObjectGroup *objGroup = NULL;
          if((rec->type == cExecObject) && 
             (rec->obj->type == cObjectGroup)) {
            objGroup = (ObjectGroup*)rec->obj;
          }
          
          switch(action) {
          case cExecutiveGroupUngroup:
            rec->group_name[0]=0;
            break;
          case cExecutiveGroupOpen:
            if(objGroup)
              objGroup->OpenOrClosed = 1;            
            break;
          case cExecutiveGroupClose:
            if(objGroup)
              objGroup->OpenOrClosed = 0;
            break;
          case cExecutiveGroupToggle:
            if(objGroup) 
              objGroup->OpenOrClosed = !objGroup->OpenOrClosed;
            break;
          case cExecutiveGroupEmpty:
            if(objGroup) {
              SpecRec *rec2 = NULL;
              while(ListIterate(I->Spec,rec2,next)) {
                if((rec2->group == rec) || WordMatchExact(G,rec2->group_name,rec->name,true)) {
                  rec2->group = NULL;
                  rec2->group_name[0] = 0;
                }
              }
            }
            break;
          case cExecutiveGroupPurge:
            if(objGroup) {
              SpecRec *rec2 = NULL;
              while(ListIterate(I->Spec,rec2,next)) {
                if((rec2->group == rec) || WordMatchExact(G,rec2->group_name,rec->name,true)) {
                  ExecutiveDelete(G,rec2->name);
                  rec2 = NULL; /* restart search (danger order N^2) */
                }
              }
            }
            break;
          case cExecutiveGroupExcise:
            if(objGroup) {

              if(rec->group_name[0]) {
                /* cascade group members up to the surrounding group */
                SpecRec *rec2 = NULL;
                while(ListIterate(I->Spec,rec2,next)) {
                  if((rec2->group == rec) ||
                     WordMatch(G,rec->name,rec2->group_name,true)) {
                    strcpy(rec2->group_name,rec->group_name);
                  }
                }
              } else if((rec->type==cExecObject)&&(rec->obj->type == cObjectGroup)) {
                /* and/or delete their group membership */
                SpecRec *rec2 = NULL;
                while(ListIterate(I->Spec,rec2,next)) {
                  if((rec2->group == rec) ||
                     WordMatch(G,rec->name,rec2->group_name,true)) {
                    rec2->group_name[0] = 0;
                  }
                }
              }
              ExecutiveDelete(G,rec->name);
            }
            break;
          }
        }
      }
      TrackerDelList(I_Tracker, list_id);
      TrackerDelIter(I_Tracker, iter_id);
      ExecutiveInvalidateGroups(G,true);
    }
  } else {
    if(obj && (obj->type==cObjectGroup)) {
      ObjectGroup *objGroup = (ObjectGroup*)obj;
      switch(action) {
      case cExecutiveGroupOpen:
        objGroup->OpenOrClosed = 1;
        break;
      case cExecutiveGroupClose:
        objGroup->OpenOrClosed = 0;
        break;
      case cExecutiveGroupToggle:
        objGroup->OpenOrClosed = !objGroup->OpenOrClosed;
        break;
      }
      if(members[0]&&(action!=cExecutiveGroupRemove))
        action = cExecutiveGroupAdd;
    
      switch(action) {
      case cExecutiveGroupAdd:
        {

          CTracker *I_Tracker= I->Tracker;
          int list_id = ExecutiveGetNamesListFromPattern(G,members,true,false);
          int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
          SpecRec *rec;
        
          while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
            if(rec && 
               ( (rec->type!=cExecObject) ||
                 ((rec->type==cExecObject) && (rec->obj!=obj)) )) {
              UtilNCopy(rec->group_name,name,sizeof(WordType));
              if(!quiet) {
                PRINTFB(G,FB_Executive, FB_Actions)
                  " Executive: adding '%s' to group '%s'.\n",rec->name,rec->group_name
                  ENDFB(G);
              }
            }
          }
          TrackerDelList(I_Tracker, list_id);
          TrackerDelIter(I_Tracker, iter_id);
        }
        break;
      }

      ExecutiveInvalidateGroups(G,true);
    }
  }
  return ok;
}

int ExecutiveGetUniqueIDObjectOffsetVLADict(PyMOLGlobals *G, 
                                            ExecutiveObjectOffset **return_vla, 
                                            OVOneToOne **return_dict)
{
  register CExecutive *I = G->Executive;
  OVOneToOne *o2o = OVOneToOne_New(G->Context->heap);
  ExecutiveObjectOffset *vla = VLAlloc(ExecutiveObjectOffset,1000);
  int n_oi = 0;
  {
    SpecRec *rec = NULL;
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type == cExecObject) {
        if(rec->obj->type == cObjectMolecule) {
          ObjectMolecule *obj = (ObjectMolecule*)rec->obj;
          register int a, id, n_atom = obj->NAtom;
          register AtomInfoType *ai = obj->AtomInfo;
          for(a=0;a<n_atom;a++) {
            if( (id=ai->unique_id) ) {
              if(OVOneToOne_GetForward(o2o,id).status == OVstatus_NOT_FOUND) {
                if(OVreturn_IS_OK(OVOneToOne_Set(o2o,id,n_oi))) {
                  VLACheck(vla,ExecutiveObjectOffset,n_oi);
                  vla[n_oi].obj = obj;
                  vla[n_oi].offset = a;
                  n_oi++;
                }
              }
            }
            ai++;
          }
        }
      }
    }
  }
  *return_dict = o2o;
  VLASize(vla,ExecutiveObjectOffset,n_oi);
  *return_vla = vla;
  return 1;
}

int ExecutiveDrawCmd(PyMOLGlobals *G, int width, int height,int antialias, int entire_window, int quiet)
{
  CExecutive *I = G->Executive;
  if((width<=0)&&(height<=0)) {
    SceneGetWidthHeight(G,&width,&height);
  }
  if(antialias<0)
    antialias = SettingGetGlobal_i(G,cSetting_antialias);
  if(entire_window) {
    SceneInvalidateCopy(G,false);
    OrthoDirty(G);
    I->CaptureFlag = true;
  } else {
    SceneDeferImage(G,width,height,NULL,antialias, -1.0, quiet);
  }
  return 1;
}

int ExecutiveMatrixCopy2(PyMOLGlobals *G,
                         CObject *source_obj, CObject *target_obj,
                         int   source_mode,  int target_mode, 
                         int   source_state, int target_state,
                         int   target_undo,
                         int   log,          int quiet)
{
  /*  mode 0: raw coordinates, as per the txf history
      mode 1: object TTT matrix
      mode 2: state matrix */

  int ok = true;
  int matrix_mode = SettingGetGlobal_b(G,cSetting_matrix_mode);
  int copy_ttt_too = false;
  if((source_mode<0)&&(target_mode<0)) {
    copy_ttt_too = true;
  }
  if(source_mode<0) 
    source_mode = matrix_mode;
  if(target_mode<0)
    target_mode = matrix_mode;
  
  switch(source_mode) {
  case 0: /* txf history is the source matrix */
    { 
      double *history = NULL;
      int found = ExecutiveGetObjectMatrix2(G,source_obj,source_state,&history,false);
      if(found) {
        switch(target_mode) {
        case 0: /* apply changes to coordinates in the target object */
          {
            double temp_inverse[16];
            if(target_undo) {
              double *target_history = NULL;
              int target_found =  ExecutiveGetObjectMatrix2(G,source_obj,
                                                            target_state,
                                                            &target_history,
                                                            false);
              if(target_found && target_history) {
                invert_special44d44d(target_history, temp_inverse);
                if(history) {
                  right_multiply44d44d(temp_inverse,history);
                  history = temp_inverse;
                } else {
                  history = temp_inverse;
                }
              }
              {
                float historyf[16];
                if(history) {
                  convert44d44f(history,historyf);
                } else {
                  identity44f(historyf);
                }
                ExecutiveTransformObjectSelection2(G,target_obj, target_state, 
                                                   "",log,historyf,true,false);
              }
            }
            if(copy_ttt_too) {
              float *tttf;
              int found = ObjectGetTTT(source_obj,&tttf,-1);
              if(found) {
                ObjectSetTTT(target_obj,tttf,-1);
                if(target_obj->fInvalidate)
                  target_obj->fInvalidate(target_obj,cRepNone,cRepInvExtents,-1);
              }
            }
          }
          break;
        case 1: /* applying changes to the object's TTT matrix */
          if(history) {
            float tttf[16];
            convertR44dTTTf(history,tttf);
            ObjectSetTTT(target_obj,tttf,-1);
          } else {
            ObjectSetTTT(target_obj,NULL,-1);
          }
          if(target_obj->fInvalidate)
            target_obj->fInvalidate(target_obj,cRepNone,cRepInvExtents,-1);
          break;
        case 2: /* applying changes to the state matrix */
          ok = ExecutiveSetObjectMatrix2(G,target_obj,target_state,history);
          break;
        }
        break;
      }
    }
    break;
  case 1: /* from the TTT matrix */
    {
      /* note that for now we're forcing states to be -1 */
      /* in the future, we may have per-state TTTs -- though right now the
         view matrices serve that purpose */
      
      float *tttf;
      int found = ObjectGetTTT(source_obj,&tttf,-1);
      if(found) {
        switch(target_mode) {
        case 0: /* coordinates & history unsupported.. */
                /* should complain */
          break;
        case 1: /* TTT */
          ObjectSetTTT(target_obj,tttf,-1);
          if(target_obj->fInvalidate)
            target_obj->fInvalidate(target_obj,cRepNone,cRepInvExtents,-1);
          break;
        case 2: /* State */
          if(tttf) {
            double homo[16];
            convertTTTfR44d(tttf,homo);
            ok = ExecutiveSetObjectMatrix2(G,target_obj,-1,homo);
          } else {
            ok = ExecutiveSetObjectMatrix2(G,target_obj,-1,NULL);
          }
          break;
        }
      }
    }
    break;
  case 2: /* from the state matrix */
    {
      double *homo;
      int found = ExecutiveGetObjectMatrix2(G,source_obj,source_state,&homo,false);
      if(found) {
        switch(target_mode) {
        case 0: /* coordinates & history */
                /* TODO */
          break;
        case 1: /* TTT */
          if(homo) {
            float tttf[16];
            convertR44dTTTf(homo,tttf);
            ObjectSetTTT(target_obj,tttf,-1);
            if(target_obj->fInvalidate)
              target_obj->fInvalidate(target_obj,cRepNone,cRepInvExtents,-1);
          } else {
            ObjectSetTTT(target_obj,NULL,-1);
            if(target_obj->fInvalidate)
              target_obj->fInvalidate(target_obj,cRepNone,cRepInvExtents,-1);
          }
          break;
        case 2: /* State */
          ok = ExecutiveSetObjectMatrix2(G,target_obj,target_state,homo);
          if(copy_ttt_too) {
            float *tttf;
            int found = ObjectGetTTT(source_obj,&tttf,-1);
            if(found) {
              ObjectSetTTT(target_obj,tttf,-1);
              if(target_obj->fInvalidate)
                target_obj->fInvalidate(target_obj,cRepNone,cRepInvExtents,-1);
            }
          }
          break;
        }
      }
    }
    break;
  }
  SceneInvalidate(G);
  return ok;
}

int ExecutiveMatrixCopy(PyMOLGlobals *G,
                        char *source_name, char *target_name,
                        int   source_mode,  int target_mode, 
                        int   source_state, int target_state,
                        int   target_undo,
                        int   log,          int quiet)
{
  /*  mode 0: raw coordinates, as per the txf history
      mode 1: object TTT matrix
      mode 2: state matrix */
  register CExecutive *I = G->Executive;
  CTracker *I_Tracker= I->Tracker;
  SpecRec *src_rec = ExecutiveFindSpec(G,source_name);
  int ok = true;
  int matrix_mode = SettingGetGlobal_b(G,cSetting_matrix_mode);
  int copy_ttt_too = false;
  if((source_mode<0)&&(target_mode<0)) {
    copy_ttt_too = true;
  }
  if(source_mode<0) 
    source_mode = matrix_mode;
  if(target_mode<0)
    target_mode = matrix_mode;
  
  switch(source_mode) {
  case 0: /* txf history is the source matrix */
    { 
      double *history = NULL;
      int found = ExecutiveGetObjectMatrix(G,source_name,source_state,&history,false);
      if(found) {

        int list_id = ExecutiveGetNamesListFromPattern(G,target_name,
                                                       true,cExecExpandKeepGroups);
        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec *rec;
        
        while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
          if(rec && (rec!=src_rec)) {
            switch(rec->type) {
            case cExecObject:

              switch(target_mode) {
              case 0: /* apply changes to coordinates in the target object */
                {
                  double temp_inverse[16];
                  if(target_undo) {
                    double *target_history = NULL;
                    int target_found =  ExecutiveGetObjectMatrix(G,rec->name,
                                                                 target_state,
                                                                 &target_history,
                                                                 false);
                    if(target_found && target_history) {
                      invert_special44d44d(target_history, temp_inverse);
                      if(history) {
                        right_multiply44d44d(temp_inverse,history);
                        history = temp_inverse;
                      } else {
                        history = temp_inverse;
                      }
                    }
                  }
                  {
                    float historyf[16];
                    if(history) {
                      convert44d44f(history,historyf);
                    } else {
                      identity44f(historyf);
                    }
                    ExecutiveTransformObjectSelection(G,rec->name, target_state, 
                                                      "",log,historyf,true,false);
                  }
                  if(copy_ttt_too) {
                    float *tttf;
                    int found = ExecutiveGetObjectTTT(G,source_name,&tttf,-1,quiet);
                    if(found) {
                      ExecutiveSetObjectTTT(G,rec->name,tttf,-1,quiet);
                    }
                  }
                }
                break;
              case 1: /* applying changes to the object's TTT matrix */
                if(history) {
                  float tttf[16];
                  convertR44dTTTf(history,tttf);
                  ExecutiveSetObjectTTT(G,rec->name,tttf,-1,quiet);
                } else {
                  ExecutiveSetObjectTTT(G,rec->name,NULL,-1,quiet);
                }
                /* to do: logging, return values, etc. */
                break;
              case 2: /* applying changes to the state matrix */
                ok = ExecutiveSetObjectMatrix(G,rec->name,target_state,history);
                break;
              }
              break;
            }
          }
        }
        TrackerDelList(I_Tracker, list_id);
        TrackerDelIter(I_Tracker, iter_id);
      }
    }
    break;
  case 1: /* from the TTT matrix */
    {
      /* note that for now we're forcing states to be -1 */
      /* in the future, we may have per-state TTTs -- though right now the
         view matrices serve that purpose */
    
      float *tttf;
      int found = ExecutiveGetObjectTTT(G,source_name,&tttf,-1,quiet);
      if(found) {

        int list_id = ExecutiveGetNamesListFromPattern(G,target_name,true,
                                                       cExecExpandKeepGroups);
        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec *rec;
        
        while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
          if(rec && (rec!=src_rec)) {

            switch(rec->type) {
            case cExecObject:
              
              switch(target_mode) {
              case 0: /* coordinates & history unsupported.. */
                /* should complain */
                break;
              case 1: /* TTT */
                ExecutiveSetObjectTTT(G,rec->name,tttf,-1,quiet);
                break;
              case 2: /* State */
                if(tttf) {
                  double homo[16];
                  convertTTTfR44d(tttf,homo);
                  ok = ExecutiveSetObjectMatrix(G,rec->name,-1,homo);
                } else {
                  ok = ExecutiveSetObjectMatrix(G,rec->name,-1,NULL);
                }
                break;
              }
            }
          }
        }

        TrackerDelList(I_Tracker, list_id);
        TrackerDelIter(I_Tracker, iter_id);
      }
    }
    break;
  case 2: /* from the state matrix */
    {
      double *homo;
      int found = ExecutiveGetObjectMatrix(G,source_name,source_state,&homo,false);
      if(found) {

        int list_id = ExecutiveGetNamesListFromPattern(G,target_name,true,cExecExpandKeepGroups);
        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        SpecRec *rec;
        
        while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
          if(rec && (rec!=src_rec)) {
            switch(rec->type) {
            case cExecObject:
        
              switch(target_mode) {
              case 0: /* coordinates & history */
                /* TODO */
                break;
              case 1: /* TTT */
                if(homo) {
                  float tttf[16];
                  convertR44dTTTf(homo,tttf);
                  ExecutiveSetObjectTTT(G,rec->name,tttf,-1,quiet);
                } else {
                  ExecutiveSetObjectTTT(G,rec->name,NULL,-1,quiet);
                }
                break;
              case 2: /* State */
                ok = ExecutiveSetObjectMatrix(G,rec->name,target_state,homo);
                if(copy_ttt_too) {
                  float *tttf;
                  int found = ExecutiveGetObjectTTT(G,source_name,&tttf,-1,quiet);
                  if(found) {
                    ExecutiveSetObjectTTT(G,rec->name,tttf,-1,quiet);
                  }
                }
                break;
              }
            }
          }
        }
        TrackerDelList(I_Tracker, list_id);
        TrackerDelIter(I_Tracker, iter_id);
      }
    }
    break;
  }
  SceneInvalidate(G);
  return ok;
}

static void ExecutiveInvalidateMapDependents(PyMOLGlobals *G,char *map_name)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      switch(rec->obj->type) {
      case cObjectMesh:
        ObjectMeshInvalidateMapName((ObjectMesh*)rec->obj,map_name);
        break;
      case cObjectSurface:
        ObjectSurfaceInvalidateMapName((ObjectSurface*)rec->obj,map_name);
        break;
      }
    }
  }
  SceneInvalidate(G);
}

void ExecutiveResetMatrix(PyMOLGlobals *G,
                          char *name,
                          int   mode,
                          int   state,
                          int   log,  
                          int quiet)
{
  register CExecutive *I = G->Executive;
  CTracker *I_Tracker= I->Tracker;
  int matrix_mode = SettingGetGlobal_b(G,cSetting_matrix_mode);
  int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;
  
  if(mode<0)
    mode = matrix_mode;
  while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
    if(rec && (rec->type==cExecObject)) {
      /*  CObject *obj = ExecutiveFindObjectByName(G,name);*/
      CObject *obj = rec->obj;
      if(obj) {
        switch(obj->type) {
        case cObjectMolecule:
          switch(mode) {
          case 0: /* transformations already applied to the coordinates */
            { 
              double *history = NULL;
              int found = ExecutiveGetObjectMatrix(G,rec->name,state,&history,false);
              if(found && history) {
                double temp_inverse[16];
                float historyf[16];
                invert_special44d44d(history, temp_inverse);
                convert44d44f(temp_inverse,historyf);
                ExecutiveTransformObjectSelection(G, rec->name, state, "",
                                                  log,historyf,true,false);
              }
            }
            break;
          case 1: /* operate on the TTT display matrix */
            ObjectResetTTT(obj);
            if(obj->fInvalidate)
              obj->fInvalidate(obj,cRepNone,cRepInvExtents,-1);
            
            break;
          case 2: /* applying changes to the state matrix */
            {
              double ident[16];
              identity44d(ident);
              ExecutiveSetObjectMatrix(G,rec->name,state,ident);
            }
            break;
          }
          break;
        case cObjectMap:
          ObjectMapResetMatrix((ObjectMap*)obj,state);
          break;
        case cObjectGroup:
          ObjectGroupResetMatrix((ObjectGroup*)obj,state);
          break;
        }
      }
    }
  }
}

static double ret_mat[16]; /* UGH ..not thread-safe */

static int ExecutiveGetObjectMatrix2(PyMOLGlobals *G,CObject *obj,int state,double **matrix, int incl_ttt)
{
  /* right now, this only makes sense for molecule objects -- but in
     time all objects should have per-state matrices */

  int ok=false;
  if(state<0) {
    /* to do -- TTT only */
  } else {
    switch(obj->type) {
    case cObjectMolecule:
      ok = ObjectMoleculeGetMatrix((ObjectMolecule*)obj,state,matrix);
      break;
    case cObjectMap:
      ok = ObjectMapGetMatrix((ObjectMap*)obj,state,matrix);
      break;
    case cObjectGroup:
      ok = ObjectGroupGetMatrix((ObjectGroup*)obj,state,matrix);
      break;
    }
    
    if(ok && incl_ttt) {
      float *ttt;
      double tttd[16];
      if(ObjectGetTTT(obj,&ttt,-1)) {
        convertTTTfR44d(ttt,tttd);
        if(*matrix) {
          copy44d(*matrix,ret_mat);
        } else {
          identity44d(ret_mat);
        }
        left_multiply44d44d(tttd,ret_mat);
        *matrix = ret_mat;
      }
    }
  }
  return ok;
}

int ExecutiveGetObjectMatrix(PyMOLGlobals *G,char *name,int state,double **matrix, int incl_ttt)
{
  int ok=false;
  CObject *obj = ExecutiveFindObjectByName(G,name);
  if(obj) {
    return ExecutiveGetObjectMatrix2(G,obj,state,matrix,incl_ttt);
  }
  return ok;
}

static int ExecutiveSetObjectMatrix2(PyMOLGlobals *G,CObject *obj,int state,double *matrix)
{
  /* -1 for the TTT matrix, 0 or greater for the state matrix */

  /* right now, this only makes sense for molecule objects -- but in
     time all objects should have per-state matrices */
  int ok=false;
  if(state<0) {
      
  } else {
    switch(obj->type) {
    case cObjectMolecule:
      ok = ObjectMoleculeSetMatrix((ObjectMolecule*)obj,state,matrix);
      break;
    case cObjectMap:
      ok = ObjectMapSetMatrix((ObjectMap*)obj,state,matrix);
      break;
    case cObjectGroup:
      ok = ObjectGroupSetMatrix((ObjectGroup*)obj,state,matrix);
      break;
    }
  }
  return ok;
}
 
int ExecutiveSetObjectMatrix(PyMOLGlobals *G,char *name,int state,double *matrix)
{
  int ok=false;
  CObject *obj = ExecutiveFindObjectByName(G,name);
  if(obj) {
    return ExecutiveSetObjectMatrix2(G,obj,state,matrix);
  }
  return ok;
}

static int ExecutiveCountNames(PyMOLGlobals *G)
{
  int count=0;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next))
    count++;
  
  return(count);
}

static int ReorderOrderFn(PyMOLGlobals *G,SpecRec **rec,int l,int r)
{
  return (WordCompare(G,rec[l]->name,rec[r]->name,true)<=0);
}

int ExecutiveOrder(PyMOLGlobals *G, char *s1, int sort,int location)
{
  register CExecutive *I = G->Executive;
  CTracker *I_Tracker= I->Tracker;
  int ok=true;
  CWordList *word_list = WordListNew(G,s1);
  int n_names = ExecutiveCountNames(G);


  if(n_names) {
    SpecRec **list,**subset,**sorted;
    int *index = NULL;
    int n_sel;
    int source_row = -1;
    list = Alloc(SpecRec*,n_names);
    subset = Calloc(SpecRec*,n_names);
    sorted = Calloc(SpecRec*,n_names);
    index = Alloc(int,n_names);
    if(list&&subset) {
      /* create an array of current names */
      {
        SpecRec *rec = NULL;
        int a = 0;
        /* copy all names into array */
        while(ListIterate(I->Spec,rec,next)) {
          list[a] = rec;
          a++;
        }
        /* unlink them */
        for(a=0;a<n_names;a++) {
          list[a]->next = NULL;
        }
      } 
#if 1
      /* transfer matching names to the subset array */
      {
        int a;
        int entry;
        int min_entry = word_list->n_word;
        char *word = NULL;
        int word_iter = 0;
        while(WordListIterate(G,word_list,&word,&word_iter)) {
          int list_id = ExecutiveGetNamesListFromPattern(G,word,true,false);
          SpecRec *rec=NULL;
          entry = word_iter-1;
          for(a=n_names-1;a>0;a--) { /* skipping zeroth */
            int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
            while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
              if(rec == list[a]) { 
                if(entry<=min_entry) {
                  source_row = a; /* where will new list be inserted...*/
                  min_entry = entry;
                }
                /* ensure that each record appears only once */
                rec->next = subset[entry];
                subset[entry] = rec;
                list[a] = NULL;
              }
            }
            TrackerDelIter(I_Tracker, iter_id);
          }
          TrackerDelList(I_Tracker, list_id);
        }
        if(word_list->n_word && WordMatchExact(G,word_list->start[0],cKeywordAll,true))
          location=-1; /* set to top if "all" is first in list */
      }
                 
#else
      /* transfer matching names to the subset array */
      {
        int a;
        int entry;
        int min_entry = word_list->n_word;
        for(a=n_names-1;a>0;a--) { /* skipping zeroth */
          entry = WordListMatch(G,word_list, list[a]->name, true);
          if(entry>=0) { /* append onto the new list */
            list[a]->next = subset[entry];
            subset[entry] = list[a];
            list[a] = NULL;
            if(entry<=min_entry) {
              source_row = a; /* takes the earliest first match */
              min_entry = entry;
            }
          }
        }
        if(word_list->n_word && WordMatchExact(G,word_list->start[0],cKeywordAll,true))
          location=-1; /* set to top if "all" is first in list */
      }
#endif
      /* expand the selected entries */
      {
        SpecRec *rec,*last;
        int b;
        n_sel = 0;
        for(b=0;b<word_list->n_word;b++) {
          rec=subset[b];
          while(rec) {
            sorted[n_sel++] = rec;
            last = rec;
            rec = rec->next;
            last->next = NULL;
          }
        }
      }
      /* sort the selected entries, if requested */
      if(sort) {
        UtilCopyMem(subset,sorted,sizeof(SpecRec*)*n_sel);
        {
          int a;
          UtilSortIndexGlobals(G,n_sel,subset,index,
                               (UtilOrderFnGlobals*)ReorderOrderFn);
          for(a=0;a<n_sel;a++) {
            sorted[a] = subset[index[a]];
          }
        }
      }
      /* reassemble the list using the new order */
      {
        SpecRec *spec= NULL;
        SpecRec *last= NULL;
        int a,b;
        int flag;
        for(a=0;a<n_names;a++) {
          flag=false;
          if(sorted) { /* not yet added */
            switch(location) {
            case -1: /* top */
              if(a==1) flag=true;
              break;
            case 0: 
              if(source_row>=0) {
                if(a==source_row)
                  flag=true;
              } else if(!list[a]) 
                flag=true;
              break;
            }
          }
          if(flag) {
            for(b=0;b<n_sel;b++) {
              if(sorted[b]) {
                if(last)
                  last->next = sorted[b];
                last = sorted[b];
                if(!spec)
                  spec = last;
              }
            }
            FreeP(sorted);
          }
          if(list[a]) {
            if(last) 
              last->next = list[a];
            last = list[a];
            if(!spec)
              spec=last;
          }
        }
        if(sorted) { /* still not yet readded? */
          for(b=0;b<n_sel;b++) {
            if(sorted[b]) {
              if(last)
                last->next = sorted[b];
              last = sorted[b];
              if(!spec)
                spec = last;
            }
          }
        }
        I->Spec=spec;
        OrthoDirty(G);
        SeqChanged(G);

      }
      
      FreeP(index);
      FreeP(sorted);
      FreeP(list);
      FreeP(subset);
    }
    ExecutiveInvalidatePanelList(G);
  }
  WordListFreeP(word_list);
  return(ok);
}

ObjectMolecule **ExecutiveGetObjectMoleculeVLA(PyMOLGlobals *G,char *sele)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  ObjectMolecule **result = NULL;
  int s1=SelectorIndexByName(G,sele);  
  if(s1>=0) {
    ObjectMoleculeOpRec op;
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_GetObjects;
    op.obj1VLA=(ObjectMolecule**)VLAlloc(CObject*,10);
    op.i1 = 0;
    ExecutiveObjMolSeleOp(G,s1,&op);
    result = (ObjectMolecule**)op.obj1VLA;
    VLASize(result,ObjectMolecule*,op.i1);
  }
  return result;
#endif
}

/* #define ExecLineHeight 18 */
#define ExecClickMargin 2
#define ExecTopMargin 0
#define ExecToggleMargin 2
#define ExecLeftMargin 1
#define ExecRightMargin 0
#define ExecToggleWidth 17
#define ExecToggleSize 16
#define ExecToggleTextShift 4

typedef struct { 
  M4XAnnoType m4x;
  ObjectMolecule *obj;
} ProcPDBRec;

int ExecutiveSetDrag(PyMOLGlobals *G,char *name, int quiet)
{
  char drag_name[] = cEditorDrag;
  int set_flag = false;
  int result = true;
  if(name[0]) {
    ObjectMolecule *obj = ExecutiveFindObjectMoleculeByName(G,name);
    if(obj) {
      SelectorCreate(G,drag_name,obj->Obj.Name,obj,true,NULL); /* for indication only */
      EditorSetDrag(G,obj,-1,quiet,SceneGetState(G));
      set_flag = true;
    } else {
      SpecRec *rec = ExecutiveFindSpec(G,name);
      if(rec) {
        if(rec->type==cExecSelection) {
          SelectorCreate(G,drag_name,name,NULL,true,NULL);
          {
            int sele = SelectorIndexByName(G,drag_name);
            obj = SelectorGetSingleObjectMolecule(G,sele);
            if(obj) {
              EditorSetDrag(G,obj,sele,quiet,SceneGetState(G));            
              set_flag = true;
            } else {
              PRINTFB(G,FB_Executive,FB_Errors)
                " Drag-Error: selection spans more than one object.\n"
                ENDFB(G);
            }
          }
        } else if(rec->type==cExecObject) {
          switch(rec->obj->type) {
          case cObjectGroup:
            PRINTFB(G,FB_Executive,FB_Errors)
              " Drag-Error: cannot drag group objects yet.\n"
              ENDFB(G);
            break;

          }
        }
      }
    }
    result = set_flag;
    if(!result) {
      EditorInactivate(G);      
      PRINTFB(G,FB_Executive,FB_Errors)
        " Drag-Error: invalid or empty selection."
        ENDFB(G);
    }
  } else {
    EditorInactivate(G);
  }
  return result;
}

int ExecutivePop(PyMOLGlobals *G,char *target,char *source,int quiet)
{
  int ok = true;
  int src;
  int result = 0;

  ExecutiveDelete(G,target);
  if(ExecutiveFindObjectMoleculeByName(G,source)) {
    ok=false;
    PRINTFB(G,FB_Executive,FB_Errors)
      " Pop-Error: source selection '%s' can't be an object.\n",source
      ENDFB(G);
    
  } else {
    src = SelectorIndexByName(G,source);
    if(src<0)
      ok=false;
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        " Pop-Error: invalid source selection name '%s'\n",source
        ENDFB(G);
    } else {
      ObjectMoleculeOpRec op;
      
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_Pop;
      SelectorCreateEmpty(G,target,true);
      op.i1 = SelectorIndexByName(G,target);
      op.i2 = 1;
      op.i3 = 0;
      ExecutiveObjMolSeleOp(G,src,&op);
      result = op.i3;
    }
  }
  if(!result) ExecutiveDelete(G,target);
  if(!ok)
    return -1;
  else
    return result;
}

int ExecutiveGetActiveAlignmentSele(PyMOLGlobals *G)
{
  char *alignment = SettingGetGlobal_s(G,cSetting_seq_view_alignment);
  int align_sele = -1;
  if( alignment && alignment[0] ) { /* explicit alignment setting name */
    align_sele = SelectorIndexByName(G,alignment);
  } else { /* otherwise, use the first active alignment */
    SpecRec *rec = NULL;
    register CExecutive *I = G->Executive;
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->visible) {
        if(rec->type==cExecObject)
          if(rec->obj->type==cObjectAlignment) {
            if(rec->obj->fUpdate) /* allow object to update selection, if necessary */
              rec->obj->fUpdate(rec->obj);
            align_sele = SelectorIndexByName(G,rec->obj->Name);
            if(align_sele>=0)
              break;
          }
      }
    }
  }
  return align_sele;
}

int ExecutiveGetActiveSele(PyMOLGlobals *G)
{
  ObjectNameType name;
  if(ExecutiveGetActiveSeleName(G,name,false,false))
    return SelectorIndexByName(G,name);
  else
    return -1;

}

int ExecutiveGetActiveSeleName(PyMOLGlobals *G,char *name, int create_new,int log)
{
  /* TODO: cache/optimize to avoid table scan */

  int result=false;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection)
      if(rec->visible) {
        strcpy(name,rec->name);
        result = true;
      }
  }
  if((!result)&&create_new) {
    if(SettingGetGlobal_b(G,cSetting_auto_number_selections)) {
      int sel_num = SettingGetGlobal_i(G,cSetting_sel_counter) + 1;
      
      SettingSetGlobal_i(G,cSetting_sel_counter,sel_num);
      sprintf(name,"sel%02d",sel_num);
      SelectorCreateEmpty(G,name,-1);
      if(log) {
        if(SettingGet(G,cSetting_logging)) {
          OrthoLineType buf2;
          sprintf(buf2,"cmd.select('%s','none')\n",name);
          PLog(G,buf2,cPLog_no_flush);
        }
      }
    } else {
      sprintf(name,"sele");
      SelectorCreateEmpty(G,name,-1);
      if(log) {
        OrthoLineType buf2;
        sprintf(buf2,"cmd.select('%s','none')\n",name);
        PLog(G,buf2,cPLog_no_flush);
      }
    }
  }
  return result;
}


int ExecutiveFixChemistry(PyMOLGlobals *G,char *s1,char *s2,int invalidate,int quiet)
{
  int sele1=SelectorIndexByName(G,s1);
  int sele2=SelectorIndexByName(G,s2);
  int ok=true;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;

  if((sele1>=0)&&(sele2>=0)) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule) {
          ObjectMoleculeFixChemistry((ObjectMolecule*)rec->obj,sele1,sele2,invalidate);
        }
    }
  }
  return ok;
}

int ExecutiveSetObjectColor(PyMOLGlobals *G,char *name,char *color,int quiet)
{
  int result = false;
  int col_ind = ColorGetIndex(G,color);
  CObject *obj = NULL;
  obj = ExecutiveFindObjectByName(G,name);
  if(obj) {
    obj->Color = col_ind;
    result = true;
  }
  return(result);
}

int ExecutiveGetObjectColorIndex(PyMOLGlobals *G,char *name)
{
  int result = -1;
  CObject *obj = NULL;
  obj = ExecutiveFindObjectByName(G,name);
  if(obj) {
    result=obj->Color;
  }
  return(result);
}

int ExecutiveGetAtomVertex(PyMOLGlobals *G,char *s1,int state,int index,float *v)
{
  int ok=false;
  int sele1 = SelectorIndexByName(G,s1);

  if(sele1>=0) {
    ok=SelectorGetSingleAtomVertex(G,sele1,state,v);
  }
  return ok;
}

int ExecutiveProcessObjectName(PyMOLGlobals *G,char *proposed,char *actual)
{
  int result = false;
  UtilNCopy(actual,proposed,sizeof(ObjectNameType));
  if(SettingGetGlobal_b(G,cSetting_validate_object_names))
    ObjectMakeValidName(actual);
  if(SettingGetGlobal_b(G,cSetting_auto_rename_duplicate_objects)) {
    if(ExecutiveValidName(G,actual)) {
      ObjectNameType candidate;
      ObjectNameType counter;
      int cnt = 2;
      while(1) {
        sprintf(counter,"_%d",cnt);
        if((strlen(actual)+strlen(counter))>=sizeof(ObjectNameType)) {
          strcpy(candidate,actual);
          candidate[sizeof(ObjectNameType)-(strlen(counter)+1)] = 0;
          strcat(candidate,counter);
        } else {
          sprintf(candidate,"%s%s",actual,counter);
        }
        if(!ExecutiveValidName(G,candidate)) {
          strcpy(actual,candidate);
          result = true;
          break;
        }
        cnt++;
      }
    }
  }
  return 1;
}

int ExecutiveSetName(PyMOLGlobals *G,char *old_name, char *new_name)
{
  int ok=true;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;
  int found = false;
  ObjectNameType name;
  UtilNCopy(name,new_name,sizeof(ObjectNameType));
  ObjectMakeValidName(name);

  if(!name[0]) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "SetName-Error: blank names not allowed.\n"
      ENDFB(G);
    ok=false;
  } else if(WordMatchExact(G,name,cKeywordSame,true)||
            SelectorNameIsKeyword(G,name)) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "SetName-Error: name '%s' is a selection keyword.\n",name
      ENDFB(G);
    ok = false;
  }
  if(ok) {   
    if(!name[0]) 
      ok=false;
    else if(!WordMatchExact(G,name,old_name,true)) {

      while(ListIterate(I->Spec,rec,next)) {
        if(found)
          break;
        switch(rec->type) {
        case cExecObject:
          if(WordMatchExact(G,rec->obj->Name,old_name,true)) {
            ExecutiveDelKey(I,rec);
            ExecutiveDelete(G,name);
            ObjectSetName(rec->obj,name);
            UtilNCopy(rec->name,rec->obj->Name,WordLength);
            ExecutiveAddKey(I,rec);          
            if(rec->obj->type == cObjectMolecule) {
              /*
                SelectorDelete(G,old_name);
                ExecutiveUpdateObjectSelection(G,rec->obj);
              */
              SelectorSetName(G,name, old_name);
              SceneChanged(G);
              SeqChanged(G);
            }
            found = true;
          }
          break;
        case cExecSelection:
          if(WordMatchExact(G,rec->name,old_name,true)) {
            if(SelectorSetName(G,name, old_name)) {
              ExecutiveDelete(G,name); /* just in case */
              ExecutiveDelKey(I,rec);
              UtilNCopy(rec->name,name,WordLength);
              ExecutiveAddKey(I,rec);
              found = true;
              OrthoDirty(G);
            }
          }
          break;
        }
      }
      if(!found)
        ok=false;
      else {
        rec=NULL;
        while(ListIterate(I->Spec,rec,next)) {
          if(WordMatchExact(G,rec->group_name,old_name,true)) {
            UtilNCopy(rec->group_name,name,WordLength);
            /* may need to rename members too... */
          }
        }
        ExecutiveInvalidateGroups(G,false);
      }
    }
  }
  return ok; 
}
#if 0
void ExecutiveLoadMOL2(PyMOLGlobals *G,CObject *origObj,char *fname,
                       char *oname, int frame, int discrete,int finish,
                       OrthoLineType buf,int multiplex,int quiet,
                       int is_string, int zoom)
{
  int ok=true;
  FILE *f;
  long size;
  char *buffer=NULL,*p;
  CObject *obj;
  char new_name[WordLength] = "";
  char *next_entry = NULL;
  int repeat_flag = true;
  int n_processed = 0;

  if(is_string) {
    buffer=fname;
  } else {
    f=fopen(fname,"rb");
    
    if(!f) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExecutiveLoadMOL2-Error: Unable to open file '%s'\n",fname
        ENDFB(G);
      ok=false;
    } else
      {
        PRINTFB(G,FB_Executive,FB_Blather)
          " ExecutiveLoadMOL2: Loading from %s.\n",fname
          ENDFB(G);
        
        fseek(f,0,SEEK_END);
        size=ftell(f);
        fseek(f,0,SEEK_SET);
        
        buffer=(char*)mmalloc(size+255);
        ErrChkPtr(G,buffer);
        p=buffer;
        fseek(f,0,SEEK_SET);
        fread(p,size,1,f);
        p[size]=0;
        fclose(f);
      }
  }

  while(repeat_flag&&ok) {
    char *start_at = buffer;
    int is_repeat_pass = false;
    int eff_frame = frame;

    if(next_entry) {
      start_at = next_entry;
      is_repeat_pass = true;
    }

    PRINTFD(G,FB_CCmd) " ExecutiveLoadMOL2-DEBUG: loading...\n" ENDFD;

    repeat_flag=false;
    next_entry = NULL;
    if(!origObj) {

      new_name[0]=0;
      obj=(CObject*)ObjectMoleculeReadMOL2Str(G,(ObjectMolecule*)origObj,
                                              start_at,eff_frame,discrete,
                                              quiet,multiplex,new_name,
                                              &next_entry);
      if(obj) {
        if(next_entry) { /* NOTE: if set, assume that multiple PDBs are present in the file */
          repeat_flag=true;
        }

        /* assign the name (if necessary) */

        if(next_entry||is_repeat_pass) {
          if(new_name[0]==0) {
            sprintf(new_name,"%s_%d",oname,n_processed+1);
          }
          ObjectSetName(obj,new_name); /* from PDB */
          ExecutiveDelete(G,new_name); /* just in case */
        } else {
          ObjectSetName(obj,oname); /* from filename/parameter */
        }

        if(obj) {
          ExecutiveManageObject(G,obj,zoom,quiet);
          if(eff_frame<0)
            eff_frame = ((ObjectMolecule*)obj)->NCSet-1;
          if(n_processed>0) {
            if(!is_string) {
              sprintf(buf," CmdLoad: loaded %d objects from \"%s\".\n",n_processed+1,fname);
            } else {
              sprintf(buf," CmdLoad: loaded %d objects from string.\n",n_processed+1);
            }
          } else {
            if(!is_string)
              sprintf(buf," CmdLoad: \"%s\" loaded as \"%s\".\n",
                      fname,oname);
            else
              sprintf(buf," CmdLoad: MOL2-string loaded into object \"%s\", state %d.\n",
                      oname,eff_frame+1);
          }
        }
      }
    } else {

      ObjectMoleculeReadMOL2Str(G,(ObjectMolecule*)origObj,
                                start_at,eff_frame,discrete,
                                quiet,multiplex,new_name,
                                &next_entry);

      if(finish) {
        ExecutiveUpdateObjectSelection(G,origObj);
        ExecutiveDoZoom(G,origObj,false,zoom,quiet);
      }
      if(eff_frame<0)
        eff_frame = ((ObjectMolecule*)origObj)->NCSet-1;
      if(!is_string) 
        sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                fname,oname,eff_frame+1);
      else
        sprintf(buf," CmdLoad: content appended into object \"%s\", state %d.\n",
                oname,eff_frame+1);
      obj = origObj;
    }

    if(obj) {
      n_processed++;
    }
  }

  if((!is_string)&&buffer) {
    mfree(buffer);
  }
  
}
#endif

int ExecutiveLoad(PyMOLGlobals *G,CObject *origObj, 
                  char *content, int content_length,
                  int content_format,
                  char *object_name, 
                  int state, int zoom, 
                  int discrete, int finish, 
                  int multiplex, int quiet, char *plugin)
{
  int ok=true;
  int is_string = false;
  int is_handled_by_python = false;

  switch(content_format) {
  case cLoadTypePDBStr:
  case cLoadTypeMOLStr:
  case cLoadTypeMMDStr:
  case cLoadTypeXPLORStr:
  case cLoadTypeMOL2Str:
  case cLoadTypeCCP4Str:
  case cLoadTypeSDF2Str:
  case cLoadTypePHIStr:
    is_string = true;
    break;
  case cLoadTypeP1M:
  case cLoadTypePMO:
  case cLoadTypeXYZ:
  case cLoadTypePDB:
  case cLoadTypeMOL:
  case cLoadTypeMMD:
  case cLoadTypeMMDSeparate:
  case cLoadTypeTOP:
  case cLoadTypeTRJ:
  case cLoadTypeCRD:
  case cLoadTypeRST:
  case cLoadTypePQR:
  case cLoadTypeMOL2:
  case cLoadTypeSDF2:
  case cLoadTypeXPLORMap:
  case cLoadTypeCCP4Map:
  case cLoadTypePHIMap:
  case cLoadTypeFLDMap:
  case cLoadTypeBRIXMap:
  case cLoadTypeGRDMap:
  case cLoadTypeDXMap:
    is_string = false;
    break;
  case cLoadTypeCUBEMap:
    is_string = true; /* this is a lie: content is actually the filename */
    break;
  case cLoadTypeCGO:
    is_string = true; /* this is a lie: contet is actually an array of floats */
    break;
  case cLoadTypePSE:
  case cLoadTypeSDF1:
  case cLoadTypeChemPyModel:
  case cLoadTypeChemPyBrick:
  case cLoadTypeChemPyMap:
  case cLoadTypeCallback:
  case cLoadTypeR3D:
    /* should never get here... */

    is_handled_by_python = true;
    break;
  }
  
  if(is_handled_by_python) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "ExecutiveLoad-Error: unable to read that file type from C\n"
      ENDFB(G);
  } else {
    OrthoLineType buf = "";
    int already_handled = false;
    
    switch(content_format) {
    case cLoadTypePDB:
    case cLoadTypePDBStr:
      {

        ok = ExecutiveProcessPDBFile(G,origObj,content,object_name,
                                     state,discrete,finish,buf,NULL,
                                     quiet,is_string,multiplex,zoom);
        /* missing return status */
      }
      already_handled = true;
      break;
    }

    if(!already_handled) {

      FILE *f;
      long size = 0;
      char *buffer=NULL,*p;
      CObject *obj = NULL;
      char new_name[WordLength] = "";
      char *next_entry = NULL;
      int repeat_flag = true;
      int n_processed = 0;
      
      if(is_string) {
        buffer = content;
        size = (long)content_length;
      } else {
        f=fopen(content,"rb");
        
        if(!f) {
          PRINTFB(G,FB_Executive,FB_Errors)
            "ExecutiveLoad-Error: Unable to open file '%s'.\n",content
            ENDFB(G);
          ok=false;
        } else {
          
          PRINTFB(G,FB_Executive,FB_Blather)
            " ExecutiveLoad: Loading from %s.\n",content
            ENDFB(G);
          
          fseek(f,0,SEEK_END);
          size=ftell(f);
          fseek(f,0,SEEK_SET);
          
          buffer=(char*)mmalloc(size+255);
          ErrChkPtr(G,buffer);
          p=buffer;
          fseek(f,0,SEEK_SET);
          fread(p,size,1,f);
          p[size]=0;
          fclose(f);
        }
      }

      while(repeat_flag&&ok) {
        char *start_at = buffer;
        int is_repeat_pass = false;
        int eff_state = state;
        int is_new = false;

        if(next_entry) {
          start_at = next_entry;
          is_repeat_pass = true;
        }
        
        PRINTFD(G,FB_CCmd) " ExecutiveLoad: loading...\n" ENDFD;
        
        repeat_flag=false;
        next_entry = NULL;

        if(!origObj) /* this is a new object */
          is_new = true;
        
        new_name[0]=0;
        
        switch(content_format) {
        case cLoadTypeMOL:
        case cLoadTypeMOLStr:
        case cLoadTypeMOL2:
        case cLoadTypeMOL2Str:
        case cLoadTypeSDF2:
        case cLoadTypeSDF2Str:
        case cLoadTypeXYZ:
        case cLoadTypeXYZStr:
          obj=(CObject*)ObjectMoleculeReadStr(G,(ObjectMolecule*)origObj,
                                              start_at,content_format,
                                              eff_state,discrete,
                                              quiet,multiplex,new_name,
                                              &next_entry);
          break;
        case cLoadTypeXPLORMap:
        case cLoadTypeXPLORStr:
          obj=(CObject*)ObjectMapLoadXPLOR(G, (ObjectMap*)origObj, start_at, eff_state, false, quiet);
          break;
        case cLoadTypePHIMap:
        case cLoadTypePHIStr:
          obj=(CObject*)ObjectMapLoadPHI(G, (ObjectMap*)origObj, start_at, eff_state, true, size, quiet);
          break;
        case cLoadTypeCCP4Map:
        case cLoadTypeCCP4Str:
          obj=(CObject*)ObjectMapLoadCCP4(G, (ObjectMap*)origObj, start_at, eff_state, true, size, quiet);
          break;
        case cLoadTypeCUBEMap:
          if(plugin) {
            obj =(CObject*)PlugIOManagerLoadVol(G, (ObjectMap*)origObj, start_at, eff_state, quiet, plugin);
          }
          break;
        case cLoadTypeCGO:
          obj=(CObject*)ObjectCGOFromFloatArray(G,(ObjectCGO*)origObj, 
                                                (float*)start_at, size, eff_state, quiet);
          break;
        }

        if(obj) {
          if(next_entry) { /* if set, then we will assume multiple objects are present,
                              and thus need to give this object its own name */
            repeat_flag=true;
          }
          
          /* assign the name (if necessary) */
          
          if(next_entry||is_repeat_pass) {
            if(is_new && (new_name[0]==0)) { /* if there wasn't a name assigned */
              sprintf(new_name,"%s_%d",object_name,n_processed+1); /* assign a default name */
            }

            ObjectSetName(obj,new_name); /* from file */
            ExecutiveDelete(G,new_name); /* just in case there is a collision */

            is_new = true; /* from now on, treat this as a new object, since indeed it is */

          } else {
            ObjectSetName(obj,object_name); /* from filename/parameter */
          }
          
          if(obj) {

            if(is_new) {
              ExecutiveManageObject(G,obj,zoom,true); /* quiet=true -- suppressing output... */
            }
            if(obj->type == cObjectMolecule) {
              if(finish) {
                ExecutiveUpdateObjectSelection(G,obj);
                ExecutiveDoZoom(G,origObj,false,zoom,quiet);
              }
            }
            switch(obj->type) {
            case cObjectMolecule:
            case cObjectMap:
              if(eff_state<0)
                eff_state = ((ObjectMolecule*)obj)->NCSet-1;
              break;
            }
            if(n_processed>0) {
              if(!is_string) {
                sprintf(buf," ExecutiveLoad: loaded %d objects from \"%s\".\n",n_processed+1,content);
              } else {
                sprintf(buf," ExecutiveLoad: loaded %d objects from string.\n",n_processed+1);
              }
            } else {
              if(!is_string)
                sprintf(buf," ExecutiveLoad: \"%s\" loaded as \"%s\", through state %d.\n",
                        content,object_name, eff_state+1);
              else
                sprintf(buf," ExecutiveLoad: content loaded into object \"%s\", through state %d.\n",
                        object_name,eff_state+1);
            }
          }
          n_processed++;
        }

      }
      if((!is_string)&&buffer) {
        mfree(buffer);
      }
    }

    if(!quiet && buf[0]) {
      PRINTFB(G,FB_Executive,FB_Actions) 
        "%s",buf
        ENDFB(G);
    }
  }
  return(ok);

}


CObject *ExecutiveGetExistingCompatible(PyMOLGlobals *G,char *oname,int type)
{
  CObject *origObj = NULL;
  origObj=ExecutiveFindObjectByName(G,oname);
  /* check for existing object of right type, delete if not */
  if(origObj) {
    int new_type = -1;
    switch(type) {
    case cLoadTypeChemPyModel:
    case cLoadTypePDB:
    case cLoadTypePDBStr:
    case cLoadTypeXYZ:
    case cLoadTypeXYZStr:
    case cLoadTypeMOL:
    case cLoadTypeMOLStr:
    case cLoadTypeMMD:
    case cLoadTypeMMDSeparate:
    case cLoadTypeMMDStr:
    case cLoadTypePMO:
    case cLoadTypeTOP:
    case cLoadTypeTRJ:
    case cLoadTypeCRD:
    case cLoadTypeMOL2:
    case cLoadTypeMOL2Str:
    case cLoadTypeSDF2:
    case cLoadTypeSDF2Str:
    case cLoadTypePQR:
    case cLoadTypeXTC:
    case cLoadTypeTRR:
    case cLoadTypeGRO:
    case cLoadTypeTRJ2:
    case cLoadTypeG96:
    case cLoadTypeDCD:
      new_type = cObjectMolecule;
      break;
    case cLoadTypeChemPyBrick:
    case cLoadTypeChemPyMap:
    case cLoadTypeXPLORMap:
    case cLoadTypeXPLORStr:
    case cLoadTypeCCP4Map:
    case cLoadTypeCCP4Str:
    case cLoadTypeFLDMap:
    case cLoadTypeGRDMap:
    case cLoadTypeDXMap:
      new_type = cObjectMap;
      break;
    case cLoadTypeCallback:
      new_type = cObjectCallback;
      break;
    case cLoadTypeCGO:
      new_type = cObjectCGO;
      break;
    }
    if (new_type!=origObj->type) {
      ExecutiveDelete(G,origObj->Name);
      origObj=NULL;
    }
  }
  return origObj;
}

int ExecutiveProcessPDBFile(PyMOLGlobals *G,CObject *origObj,char *fname,
                            char *oname, int frame, int discrete,int finish,
                            OrthoLineType buf,PDBInfoRec *pdb_info,int quiet,
                            int is_string,int multiplex,int zoom)
{
  int ok=true;
  FILE *f;
  long size;
  char *buffer=NULL,*p;
  CObject *obj;
  char pdb_name[WordLength] = "";
  char cur_name[WordLength] = "";
  char *next_pdb = NULL;
  int repeat_flag = true;
  ProcPDBRec *processed= NULL;
  int n_processed = 0;
  int m4x_mode = 0; /* 0 = annotate, 1 = alignment */
  ProcPDBRec *target_rec = NULL;
  char nbrhood_sele[] = "m4x_nearby";
  ProcPDBRec *current = NULL;
  PDBInfoRec pdb_info_rec;
  int model_number;

  if(!pdb_info) {
    UtilZeroMem(&pdb_info_rec,sizeof(PDBInfoRec));
    pdb_info=&pdb_info_rec;
  }
  pdb_info->multiplex = multiplex;
  if(is_string) {
    buffer=fname;
  } else {
    f=fopen(fname,"rb");
    
    if(!f) {
      PRINTFB(G,FB_ObjectMolecule,FB_Errors)
        "ExecutiveProcessPDBFile-Error: Unable to open file '%s'.\n",fname
        ENDFB(G);
      ok=false;
    } else
      {
        PRINTFB(G,FB_ObjectMolecule,FB_Blather)
          " ExecutiveProcessPDBFile: Loading from %s.\n",fname
          ENDFB(G);
        
        fseek(f,0,SEEK_END);
        size=ftell(f);
        fseek(f,0,SEEK_SET);
        
        buffer=(char*)mmalloc(size+255);
        ErrChkPtr(G,buffer);
        p=buffer;
        fseek(f,0,SEEK_SET);
        fread(p,size,1,f);
        p[size]=0;
        fclose(f);
      }
  }

  if(ok) {
    processed = VLACalloc(ProcPDBRec,10);
  }
  while(repeat_flag&&ok) {
    char *start_at = buffer;
    int is_repeat_pass = false;
    int eff_frame = frame;
    int is_new = false;
    CObject *tmpObj;

    VLACheck(processed,ProcPDBRec,n_processed);
    current = processed + n_processed;

    PRINTFD(G,FB_CCmd) " ExecutiveProcessPDBFile-DEBUG: loading PDB\n" ENDFD;

    if(next_pdb) {
      start_at = next_pdb;
      is_repeat_pass = true;
    }

    M4XAnnoInit(&current->m4x);

    repeat_flag=false;
    next_pdb = NULL;
    if(!origObj) {

      is_new = true;
      pdb_name[0]=0;
      model_number = 0;
      obj=(CObject*)ObjectMoleculeReadPDBStr(G,(ObjectMolecule*)origObj,
                                             start_at,eff_frame,discrete,
                                             &current->m4x,pdb_name,
                                             &next_pdb,pdb_info,quiet, 
                                             &model_number);

    } else {
      model_number = 0;
      ObjectMoleculeReadPDBStr(G,(ObjectMolecule*)origObj,
                               start_at,eff_frame,discrete,&current->m4x,
                               pdb_name,&next_pdb,pdb_info,quiet, &model_number);
      
      
      if(finish) {
        ExecutiveUpdateObjectSelection(G,origObj);
        ExecutiveDoZoom(G,origObj,false,zoom,quiet);
      }
      if(eff_frame<0)
        eff_frame = ((ObjectMolecule*)origObj)->NCSet-1;
      if(buf) {
        if(!is_string) 
          sprintf(buf," CmdLoad: \"%s\" appended into object \"%s\", state %d.\n",
                  fname,oname,eff_frame+1);
        else
          sprintf(buf," CmdLoad: PDB-string appended into object \"%s\", state %d.\n",
                  oname,eff_frame+1);
      }
      obj = origObj;
    }

    if(obj) {
      if(next_pdb) { /* NOTE: if set, assume that multiple PDBs are present in the file */
        repeat_flag=true;
      }
    }

    if(is_new) {
      if(obj) {
        if(current->m4x.xname_flag) { /* USER XNAME trumps the PDB Header name */                 
          ObjectSetName(obj,current->m4x.xname); /* from PDB */
          if( (tmpObj = ExecutiveFindObjectByName(G,obj->Name))) {
            if(tmpObj->type != cObjectMolecule)
              ExecutiveDelete(G,current->m4x.xname); /* just in case */
            else {
              if(is_repeat_pass) {  /* this is a workaround for when PLANET accidentally duplicates the target */
                {
                  int a;
                  for(a=0; a<n_processed; a++) {
                    ProcPDBRec *cur = processed + a;
                    if(cur->obj == (ObjectMolecule*)tmpObj) {
                      cur->m4x.invisible = false;
                    }
                  }
                }
                ObjectMoleculeFree((ObjectMolecule*)obj);
                obj=NULL;
              }
            }
          }
        } else if(next_pdb) {
          if(pdb_name[0]==0) {
            if(cur_name[0]) {
              sprintf(pdb_name,"%s_%04d",cur_name,n_processed+1);
            } else {
              sprintf(pdb_name,"%s_%04d",oname,n_processed+1);
            }
          } else if(multiplex>0) {
            if(pdb_info->multi_object_status == 1 ) { /* this is a multi-object PDB file */
              strcpy(cur_name,pdb_name);
            } else if(cur_name[0]==0) {
              strcpy(cur_name, oname);
            }
            if(model_number>0) {
              sprintf(pdb_name,"%s_%04d",cur_name,model_number);
            } else {
              sprintf(pdb_name,"%s_%04d",cur_name,n_processed+1);
            }
          }
          ObjectSetName(obj,pdb_name); 
          ExecutiveDelete(G,pdb_name); /* just in case */
        } else {
          if(is_repeat_pass) {
            if(pdb_name[0]==0) {
              if(cur_name[0]) {
                sprintf(pdb_name,"%s_%04d",cur_name,n_processed+1);
              } else {
                sprintf(pdb_name,"%s_%04d",oname,n_processed+1);
              }
            } else if(multiplex>0) { 
              if(pdb_info->multi_object_status == 1 ) { /* this is a multi-object PDB file */
                strcpy(cur_name,pdb_name);
              } else if(cur_name[0]==0) {
                strcpy(cur_name, oname);
              }
              if(model_number>0) {
                sprintf(pdb_name,"%s_%04d",cur_name,model_number);
              } else {
                sprintf(pdb_name,"%s_%04d",cur_name,n_processed+1);
              }
            }
            ObjectSetName(obj,pdb_name); /* from PDB */
            ExecutiveDelete(G,pdb_name); /* just in case */
          } else {
            ObjectSetName(obj,oname); /* from filename/parameter */
          }
        }

        if(obj) {
          ExecutiveManageObject(G,obj,zoom,true);
          if(eff_frame<0)
            eff_frame = ((ObjectMolecule*)obj)->NCSet-1;
          if(buf) {
            if(n_processed<1) {
              if(!is_string)
                sprintf(buf," CmdLoad: \"%s\" loaded as \"%s\".\n",
                        fname,oname);
              else
                sprintf(buf," CmdLoad: PDB-string loaded into object \"%s\", state %d.\n",
                        oname,eff_frame+1);
            } else {
              if(!is_string) {
                sprintf(buf," CmdLoad: loaded %d objects from \"%s\".\n",n_processed+1,fname);
              } else {
                sprintf(buf," CmdLoad: loaded %d objects from string.\n",n_processed+1);
              }
            }
          }
            
        }
      }
    } 

    if(obj&&current) {
      current->obj = (ObjectMolecule*)obj;
      n_processed++;
    }
  }

  /* BEGIN METAPHORICS ANNOTATION AND ALIGNMENT CODE */

  /* sanity check -- make sure all objects are present */
  if(ok&&n_processed) {
    int a;
    for(a=0; a<n_processed; a++) {
      ProcPDBRec *current = processed + a;
      if(!ExecutiveValidateObjectPtr(G,(CObject*)current->obj,cObjectMolecule)) {
        PRINTFB(G,FB_Executive,FB_Errors)
          " Error: Missing object! possible invalid/corrupt file.\n"
          ENDFB(G);
        ok=false;
        break;
      }
    }
  }
  
  if(ok&&n_processed) { /* first, perform any Metaphorics alignment */
    /* is there a target structure? */
    {
      int a;
      for(a=0; a<n_processed; a++) {
        ProcPDBRec *current = processed + a;
        M4XAnnoType *m4x = &current->m4x;

        if(m4x->annotated_flag && m4x->align) {
          if(WordMatchExact(G,current->obj->Obj.Name,m4x->align->target,true)) {
            target_rec = current;
            break;
          }
        }
      }
    }
    if(target_rec) { /* there is a target.. */
      
      /* first, convert all IDs to genuine atom indices */

      {
        int a;
        for(a=0; a<n_processed; a++) {
          ProcPDBRec *current = processed + a;
          M4XAnnoType *m4x = &current->m4x;
          if(m4x->align) {
            ObjectMoleculeConvertIDsToIndices(current->obj, 
                                              m4x->align->id_at_point,
                                              m4x->align->n_point);
          }
        }
      }

      /* next, peform the alignments against the target */

      {
        int a;
        char aligned_name[] = "m4x_aligned";
        char tmp_sele[WordLength*4];

        SelectorCreateEmpty(G,"m4x_aligned",-1);
                
                
        for(a=0; a<n_processed; a++) {
          ProcPDBRec *current = processed + a;
          if( current != target_rec ) {
            M4XAnnoType *m4x = &current->m4x;
            if(m4x->align) {
              ObjMolPairwise pairwise;

              m4x_mode = 1; /* performing an alignment */

              ObjMolPairwiseInit(&pairwise);
              pairwise.trg_obj = target_rec->obj;
              pairwise.mbl_obj = current->obj;
              
              /* create ordered matches */

              {
                M4XAlignType *trg_align = target_rec->m4x.align;
                M4XAlignType *mbl_align = m4x->align;

                int n_point;
                int a;
                  
                n_point = trg_align->n_point;
                if(n_point > mbl_align->n_point) n_point = mbl_align->n_point;
                
                VLACheck(pairwise.trg_vla,int,n_point);
                VLACheck(pairwise.mbl_vla,int,n_point);
                
                for(a=0;a<n_point;a++) {
                  int trg_index = trg_align->id_at_point[a];
                  int mbl_index = mbl_align->id_at_point[a];
                  if((trg_index>=0)&&(mbl_index>=0)&&
                     (mbl_align->fitness[a]>=0.0F)) {
                    pairwise.trg_vla[pairwise.n_pair] = trg_index;
                    pairwise.mbl_vla[pairwise.n_pair] = mbl_index;
                    pairwise.n_pair++;
                  }
                }

                {
                  char trg_sele[WordLength],mbl_sele[WordLength];
                  char align_name[WordLength];
                  SelectorGetUniqueTmpName(G,trg_sele);
                  SelectorGetUniqueTmpName(G,mbl_sele);
                  
                  SelectorCreateOrderedFromObjectIndices(G,trg_sele,pairwise.trg_obj,
                                                         pairwise.trg_vla,pairwise.n_pair);
                  SelectorCreateOrderedFromObjectIndices(G,mbl_sele,pairwise.mbl_obj,
                                                         pairwise.mbl_vla,pairwise.n_pair);
                  
                  sprintf(align_name,"%s_%s_alignment",
                          pairwise.trg_obj->Obj.Name,
                          pairwise.mbl_obj->Obj.Name);

                  ExecutiveRMS(G,mbl_sele, trg_sele, 2, 0.0F, 0, 0,
                               align_name, 0, 0, true, 0, NULL);
                  ExecutiveColor(G,align_name,"white",0,true);
                  if(target_rec->m4x.invisible) 
                    sprintf(tmp_sele, "(%s) | (%s)",aligned_name,mbl_sele);
                  else 
                    sprintf(tmp_sele, "(%s) | (%s) | (%s)",aligned_name, trg_sele, mbl_sele);
                  SelectorCreateSimple(G,aligned_name, tmp_sele);
                  /*ExecutiveDelete(G,trg_sele);
                    ExecutiveDelete(G,mbl_sele);*/
                }

              }
              ObjMolPairwisePurge(&pairwise);
            }
          }
        }
        sprintf(tmp_sele, "bychain %s", aligned_name);
        SelectorCreateSimple(G,aligned_name, tmp_sele);
        sprintf(tmp_sele, "byres (%s expand 3.5)", aligned_name);
        SelectorCreateSimple(G,nbrhood_sele, tmp_sele);
      }
    }
  }
  
  if(ok&&n_processed) { /* next, perform any and all Metaphorics annotations */
    int a;
    int nbr_sele = SelectorIndexByName(G,nbrhood_sele);
      
    for(a=0; a<n_processed; a++) {
      ProcPDBRec *current = processed + a;
      if(current->m4x.annotated_flag) {
        char annotate_script[] = "@$PYMOL_SCRIPTS/metaphorics/annotate.pml";
        char align_script[] = "@$PYMOL_SCRIPTS/metaphorics/alignment.pml";
        char *script_file = NULL;
          
        if(a==(n_processed-1)) { /* for multi-PDB files, don't execute script until after the last file */
          switch(m4x_mode) {
          case 0:
            script_file = annotate_script;
            break;
          case 1:
            script_file = align_script;
            break;
          }
        }

        if((current!=target_rec)||(!current->m4x.invisible)) /* suppress annotations if target invisible */
          ObjectMoleculeM4XAnnotate(current->obj,&current->m4x,script_file,(m4x_mode==1),
                                    nbr_sele);
        M4XAnnoPurge(&current->m4x);
      }
    }
  }
  /* END METAPHORICS ANNOTATION AND ALIGNMENT CODE */
  
  VLAFreeP(processed);
  if((!is_string)&&buffer) {
    mfree(buffer);
  }
 
  return ok;
}

int  ExecutiveAssignSS(PyMOLGlobals *G,char *target,int state,char *context,int preserve,int quiet)
{
  int sele0=-1;
  int sele1=-1;
  int ok = false;
  sele0 = SelectorIndexByName(G,target);
  if(sele0>=0) {
    if((!context)||(!context[0])) {
      sele1=sele0;
    } else {
      sele1 = SelectorIndexByName(G,context);
    }
    if(sele1>=0) {
      ok =  SelectorAssignSS(G,sele0,sele1,state,preserve,quiet);
    }
  }
  return(ok);
}

PyObject *ExecutiveGetVisAsPyDict(PyMOLGlobals *G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result=NULL,*list,*repList;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int a;
  int n_vis;
  result = PyDict_New();
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->name[0]!='_') {
      list = PyList_New(4);
      PyList_SetItem(list,0,PyInt_FromLong(rec->visible));

      /* all executive entries have repOn */
      n_vis=0;
      for(a=0;a<cRepCnt;a++) {
        if(rec->repOn[a])
          n_vis++;
      }
      repList = PyList_New(n_vis);
      n_vis=0;
      for(a=0;a<cRepCnt;a++) {
        if(rec->repOn[a]) {
          PyList_SetItem(repList,n_vis,PyInt_FromLong(a));
          n_vis++;
        }
      }
      PyList_SetItem(list,1,repList);
      
      if(rec->type!=cExecObject) {
        PyList_SetItem(list,2,PConvAutoNone(Py_None));
        PyList_SetItem(list,3,PConvAutoNone(Py_None));
      } else { 
        /* objects have their own visib list too */
        n_vis=0;
        for(a=0;a<cRepCnt;a++) {
          if(rec->obj->RepVis[a])
            n_vis++;
        }
        repList = PyList_New(n_vis);
        n_vis=0;
        for(a=0;a<cRepCnt;a++) {
          if(rec->obj->RepVis[a]) {
            PyList_SetItem(repList,n_vis,PyInt_FromLong(a));
            n_vis++;
          }
        }
        PyList_SetItem(list,2,repList);
        PyList_SetItem(list,3,PyInt_FromLong(rec->obj->Color));
      }

      PyDict_SetItemString(result,rec->name,list);
      Py_DECREF(list);
    }
  }
  return(result);
#endif
}

int ExecutiveSetVisFromPyDict(PyMOLGlobals *G,PyObject *dict)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  WordType name;
  PyObject *key,*list,*col;
  PyObject *vis_list = NULL;
#if (PY_MAJOR_VERSION>=2)&&(PY_MINOR_VERSION>=5)
  Py_ssize_t pos = 0;
#else
  int pos = 0;
#endif
  SpecRec *rec;
  int n_vis;
  int rep;
  int a;
  int ll=0;
  if(ok) ok=(dict!=NULL);
  if(ok) ok=PyDict_Check(dict);
  if(ok) {

    SceneObjectDel(G,NULL); /* remove all objects from scene */
    ExecutiveInvalidateSceneMembers(G);

    while (PyDict_Next(dict, &pos, &key, &list)) {
      if(!PConvPyStrToStr(key,name,sizeof(WordType))) {
        ok=false;
      } else {
        rec = ExecutiveFindSpec(G,name);
        if(rec) {
          if(ok) ok = (list!=NULL);
          if(ok) ok = PyList_Check(list);
          if(ok) ll = PyList_Size(list);
          if(ok) ok = (ll>=2);
          if(ok) ok = PConvPyObjectToInt(PyList_GetItem(list,0),&rec->visible);
          if(ok) { /* rec visibility */
            vis_list = PyList_GetItem(list,1);
            if(ok) ok = (vis_list!=NULL);
            if(ok) ok = PyList_Check(vis_list);
            if(ok) {
              n_vis = PyList_Size(vis_list);
              for(a=0;a<cRepCnt;a++)
                rec->repOn[a]=false;
              for(a=0;a<n_vis;a++) {
                if(PConvPyObjectToInt(PyList_GetItem(vis_list,a),&rep)) {
                  if((rep>=0)&&(rep<cRepCnt))
                    rec->repOn[rep]=true;
                }
              }
            }
          }

          if(ok&&(rec->type==cExecObject)) { /* object properties */

            if(ll>2) { /* object visibility */
              vis_list = PyList_GetItem(list,2);
              if(ok) ok = (vis_list!=NULL);
              if(ok) if(PyList_Check(vis_list)) {
                n_vis = PyList_Size(vis_list);
                for(a=0;a<cRepCnt;a++)
                  rec->obj->RepVis[a]=false;
                for(a=0;a<n_vis;a++) {
                  if(PConvPyObjectToInt(PyList_GetItem(vis_list,a),&rep)) {
                    if((rep>=0)&&(rep<cRepCnt))
                      rec->obj->RepVis[rep]=true;
                  }
                }
              }
            }
            if(ll>3) { /* object color */
              col = PyList_GetItem(list,3);
              if(ok) ok = (col!=NULL);
              if(ok) if(PyInt_Check(col)) {
                ok = PConvPyObjectToInt(col,&rec->obj->Color);
                if(rec->obj->fInvalidate)
                  rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvColor,-1);
              }
            }
          }
          if(rec->visible&&(rec->type==cExecObject)) {
            rec->in_scene = SceneObjectAdd(G,rec->obj); 
            ExecutiveInvalidateSceneMembers(G);
          }
        }
      }
    }
  }
  return ok;
#endif
}

int ExecutiveIsolevel(PyMOLGlobals *G,char *name,float level,int state,int query,float *result,int quiet)
{
  int ok =true;
  CObject *obj;
  obj = ExecutiveFindObjectByName(G,name);
  if(obj) {
    switch(obj->type) {
    case cObjectMesh:
      if(!query) {
        ObjectMeshSetLevel((ObjectMesh*)obj,level,state,quiet);
        SceneChanged(G);
      } else if(result) {
        ok = ObjectMeshGetLevel((ObjectMesh*)obj,state,result);
      }
      break;
    case cObjectSurface:
      if(!query) {
        ObjectSurfaceSetLevel((ObjectSurface*)obj,level,state,quiet);
        SceneChanged(G);
      } else if(result) {
        ok = ObjectSurfaceGetLevel((ObjectSurface*)obj,state,result);
      }
      break;
    default:
      ok=false;
      PRINTFB(G,FB_Executive,FB_Errors)
        " Isolevel-Error: object \"%s\" is of wrong type.",name
        ENDFB(G);
      break;
    }
  }
  return(ok);

}

int ExecutiveSpectrum(PyMOLGlobals *G,char *s1,char *expr,float min,float max,int first,int last,
                      char *prefix,int digits,int byres,int quiet,
                      float *min_ret,float *max_ret)
{
  int ok=true;
  int sele1;
  int n_color,n_atom;
  ObjectMoleculeOpRec op;
  WordType buffer;
  int *color_index = NULL;
  float *value = NULL;
  int a,b;
  char pat[] = "%0Xd";
  int pref_len;
  char *at;
  float range;

  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {

    if(digits>9) digits = 9;
    pat[2]=('0'+digits);
    UtilNCopy(buffer,prefix,sizeof(WordType)-digits);
    
    pref_len = strlen(prefix);
    at = buffer+pref_len;

    n_color = abs(first-last)+1;
    if(n_color) {
      color_index = Alloc(int,n_color);
      for(a=0;a<n_color;a++) {
        b = first + ((last-first)*a)/(n_color-1);
        sprintf(at,pat,b);
        color_index[a] = ColorGetIndex(G,buffer);
      }
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_CountAtoms;
      op.i1=0;
      ExecutiveObjMolSeleOp(G,sele1,&op);
      n_atom=op.i1;
      
      if(n_atom) {
        value = Calloc(float,n_atom);
        
        if(WordMatch(G,"count",expr,true)) {
          for(a=0;a<n_atom;a++) {
            value[a]=(float)a+1;
          }
        } else if(WordMatch(G,"b",expr,true)) {
          op.code = OMOP_GetBFactors;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(G,sele1,&op);
        } else if(WordMatch(G,"q",expr,true)) {
          op.code = OMOP_GetOccupancies;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(G,sele1,&op);
        } else if(WordMatch(G,"pc",expr,true)) {
          op.code = OMOP_GetPartialCharges;
          op.i1 = 0;
          op.ff1 = value;
          ExecutiveObjMolSeleOp(G,sele1,&op);
        }

        if(max<min) {
          max = value[0];
          min = value[0];
          for(a=1;a<n_atom;a++) {
            if(value[a]<min) min=value[a];
            if(value[a]>max) max=value[a];
          }
        }
        range = max-min;

        if(!quiet) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Spectrum: range (%8.5f to %8.5f).\n"
            ,min,max
            ENDFB(G);
        }
        if(range==0.0F)
          range = 1.0F;
        *min_ret = min;
        *max_ret = max;


        op.code = OMOP_Spectrum;
        op.i1 = n_color-1;
        op.i2 = n_atom;
        op.i3 = 0;
        op.i4 = byres;
        op.ii1 = color_index;
        op.ff1 = value;
        op.f1 = min;
        op.f2 = range;
        
        ExecutiveObjMolSeleOp(G,sele1,&op);

        op.code=OMOP_INVA;
        op.i1=cRepAll; 
        op.i2=cRepInvColor;
        ExecutiveObjMolSeleOp(G,sele1,&op);

      }
    }

    FreeP(color_index);
    FreeP(value);
  }
  return(ok);
}

char *ExecutiveGetChains(PyMOLGlobals *G,char *sele,int state,int *null_chain)
{
  int sele1;
  char *result = NULL;
  int chains[256];
  int a,c;
  ObjectMoleculeOpRec op;

  sele1 = SelectorIndexByName(G,sele);
  if(sele1>=0) {

    for(a=0;a<256;a++) {
      chains[a]=0;
    }
    ObjectMoleculeOpRecInit(&op);
    op.code=OMOP_GetChains;
    op.ii1 = chains;
    op.i1=0;
    ExecutiveObjMolSeleOp(G,sele1,&op);
    c=0;
    for(a=1;a<256;a++) {
      if(chains[a]) c++;
    }
    result = Calloc(char,c+1);
    if(result) {
      c=0;
      *null_chain = chains[0];
      for(a=1;a<256;a++) {
        if(chains[a]) {
          result[c]=(char)a;
          c++;
        }
      }
    }
  } else {
    ErrMessage(G,"ExecutiveGetChains","Bad selection.");
  }
  return(result);
}

int ExecutiveValidateObjectPtr(PyMOLGlobals *G,CObject *ptr,int object_type)
{ 
  /* this routine needs to be sped up significantly... */

  register CExecutive *I = G->Executive;
  int ok=false;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->obj == ptr) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==object_type) {
          ok=true;
          break;
        }
      }
    }
  }
  return(ok);
}

int ExecutiveRampNew(PyMOLGlobals *G,char *name,char *src_name,
                     float *range,float *color,
                     int src_state,char *sele, float beyond,
                     float within,float sigma,int zero,int calc_mode,int quiet)
{
  ObjectGadgetRamp *obj = NULL;
  int ok =true;
  CObject *src_obj = NULL;
  float *vert_vla = NULL;
  src_obj = ExecutiveFindObjectByName(G,src_name);
  if(src_obj) {
    if((src_obj->type!=cObjectMap)&&
       (src_obj->type!=cObjectMolecule)) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExecutiveRampMapNew: Error: object '%s' is not a map or molecule.\n",src_name
        ENDFB(G);
      ok=false;
    }
  } else if(WordMatch(G,src_name,cKeywordNone,true)) {
    src_obj=NULL;
  } else {
    PRINTFB(G,FB_Executive,FB_Errors)
      "ExecutiveRampMapNew: Error: object '%s' not found.\n",src_name
      ENDFB(G);
    ok = false;
  }
  if(ok) {
    if(!src_obj) {
      ok = ok && (obj=ObjectGadgetRampMolNewAsDefined(G,NULL,range,color,src_state,calc_mode));
    } else {
      switch(src_obj->type) {
      case cObjectMap:
        if(sele&&sele[0]) {
          vert_vla = ExecutiveGetVertexVLA(G,sele,src_state);
        }
        ok = ok && (obj=ObjectGadgetRampMapNewAsDefined(G,(ObjectMap*)src_obj,
                                                        range,color,src_state,
                                                        vert_vla,beyond,within,
                                                        sigma,zero,calc_mode));
        break;
      case cObjectMolecule:
        ok = ok && (obj=ObjectGadgetRampMolNewAsDefined(G,(ObjectMolecule*)src_obj,
                                                        range,color,src_state,calc_mode));
        break;
      }
    }
  }
  if(ok) ExecutiveDelete(G,name); 
  if(ok) ObjectSetName((CObject*)obj,name);
  if(ok) ColorRegisterExt(G,name,(void*)obj,cColorGadgetRamp);
  if(ok) ExecutiveManageObject(G,(CObject*)obj,false,quiet);
  if(ok) ExecutiveInvalidateRep(G,cKeywordAll,cRepAll,cRepInvColor); /* recolor everything */
  VLAFreeP(vert_vla);
  return(ok);
}


#ifndef _PYMOL_NOPY
static PyObject *ExecutiveGetExecObjectAsPyList(PyMOLGlobals *G,SpecRec *rec)
{

  PyObject *result = NULL;
  result = PyList_New(7);
  PyList_SetItem(result,0,PyString_FromString(rec->obj->Name));
  PyList_SetItem(result,1,PyInt_FromLong(cExecObject));
  PyList_SetItem(result,2,PyInt_FromLong(rec->visible));
  PyList_SetItem(result,3,PConvIntArrayToPyList(rec->repOn,cRepCnt));
  PyList_SetItem(result,4,PyInt_FromLong(rec->obj->type));
  switch(rec->obj->type) {
  case cObjectGadget:
    PyList_SetItem(result,5,ObjectGadgetAsPyList((ObjectGadget*)rec->obj));    
    break;
  case cObjectMolecule:
    PyList_SetItem(result,5,ObjectMoleculeAsPyList((ObjectMolecule*)rec->obj));
    break;
  case cObjectMeasurement:
    PyList_SetItem(result,5,ObjectDistAsPyList((ObjectDist*)rec->obj));
    break;
  case cObjectMap:
    PyList_SetItem(result,5,ObjectMapAsPyList((ObjectMap*)rec->obj));
    break;
  case cObjectMesh:
    PyList_SetItem(result,5,ObjectMeshAsPyList((ObjectMesh*)rec->obj));
    break;
  case cObjectSlice:
    PyList_SetItem(result,5,ObjectSliceAsPyList((ObjectSlice*)rec->obj));
    break;
  case cObjectSurface:
    PyList_SetItem(result,5,ObjectSurfaceAsPyList((ObjectSurface*)rec->obj));
    break;
  case cObjectCGO:
    PyList_SetItem(result,5,ObjectCGOAsPyList((ObjectCGO*)rec->obj));
    break;
  case cObjectAlignment:
    PyList_SetItem(result,5,ObjectAlignmentAsPyList((ObjectAlignment*)rec->obj));
    break;
  case cObjectGroup:
    PyList_SetItem(result,5,ObjectGroupAsPyList((ObjectGroup*)rec->obj));
    break;
  default: 
    PyList_SetItem(result,5,PConvAutoNone(NULL));
    break;
  }
  PyList_SetItem(result,6,PyString_FromString(rec->group_name)); 

  return(result);  
}

static int ExecutiveSetNamedEntries(PyMOLGlobals *G,PyObject *names,int version,
                                    int part_rest,int part_sess)
{
  register CExecutive *I = G->Executive;  
  int ok=true;
  int skip=false;
  int a=0,l=0,ll=0;
  PyObject *cur;
  SpecRec *rec = NULL;
  int extra_int;
  int incomplete = false;

  if(ok) ok = (names!=NULL);
  if(ok) ok = PyList_Check(names);
  if(ok) l = PyList_Size(names);
  while(ok&&(a<l)) {
    cur = PyList_GetItem(names,a);
    if(cur!=Py_None) { /* skip over None w/o aborting */
      skip=false;
      rec=NULL;
      ListElemCalloc(G,rec,SpecRec); 
      rec->next=NULL;
      rec->name[0]=0;
      if(ok) ok = PyList_Check(cur);
      if(ok) ll = PyList_Size(cur);
      if(ok) ok = PConvPyStrToStr(PyList_GetItem(cur,0),rec->name,sizeof(WordType));
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,1),&rec->type);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,2),&rec->visible);
      if(ok) ok = PConvPyListToIntArrayInPlaceAutoZero(PyList_GetItem(cur,3),
                                                       rec->repOn,cRepCnt);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,4),&extra_int); 
      switch(rec->type) {
      case cExecObject:
        switch(extra_int) {
        case cObjectMolecule:
          if(ok) ok = ObjectMoleculeNewFromPyList(G,PyList_GetItem(cur,5),(ObjectMolecule**)&rec->obj);
          break;
        case cObjectMeasurement:
          if(ok) ok = ObjectDistNewFromPyList(G,PyList_GetItem(cur,5),(ObjectDist**)&rec->obj);
          break;
        case cObjectMap:
          if(ok) ok = ObjectMapNewFromPyList(G,PyList_GetItem(cur,5),(ObjectMap**)&rec->obj);
          break;
        case cObjectMesh:
          if(ok) ok = ObjectMeshNewFromPyList(G,PyList_GetItem(cur,5),(ObjectMesh**)&rec->obj);
          break;
        case cObjectSlice:
          if(ok) ok = ObjectSliceNewFromPyList(G,PyList_GetItem(cur,5),(ObjectSlice**)&rec->obj);
          break;	  
        case cObjectSurface:
          if(ok) ok = ObjectSurfaceNewFromPyList(G,PyList_GetItem(cur,5),(ObjectSurface**)&rec->obj);
          break;
        case cObjectCGO:
          if(ok) ok = ObjectCGONewFromPyList(G,PyList_GetItem(cur,5),(ObjectCGO**)&rec->obj,version);
          break;
        case cObjectGadget:
          if(ok) ok = ObjectGadgetNewFromPyList(G,PyList_GetItem(cur,5),(ObjectGadget**)&rec->obj,version);
          break;
        case cObjectAlignment:
          if(ok) ok = ObjectAlignmentNewFromPyList(G,PyList_GetItem(cur,5),(ObjectAlignment**)&rec->obj,version);
          break;
        case cObjectGroup:
          if(ok) ok = ObjectGroupNewFromPyList(G,PyList_GetItem(cur,5),(ObjectGroup**)&rec->obj,version);
          break;
          
        default:
          PRINTFB(G,FB_Executive,FB_Errors)
            " Executive: skipping unrecognized object \"%s\" of type %d.\n",
            rec->name,extra_int
            ENDFB(G);
          skip=true;
          break;
        }
        break;
      case cExecSelection: /* on the first pass, just create an entry in the rec list */
        rec->sele_color=extra_int;
        if(part_rest||part_sess) { /* don't attempt to restore selections with partial sessions */
          skip=true;
        }
        break;
      }

      if(ll>6) {
        if(ok) ok = PConvPyStrToStr(PyList_GetItem(cur,6),rec->group_name,sizeof(WordType));    
      }

      if(PyErr_Occurred()) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetNamedEntries-Error: after object \"%s\".\n",rec->name
          ENDFB(G);
        PyErr_Print();  
      }

      if(ok&&!skip) {
        switch(rec->type) {
        case cExecObject:        
          if(rec->visible) {
            rec->in_scene = SceneObjectAdd(G,rec->obj);
            ExecutiveInvalidateSceneMembers(G);
          }
          ExecutiveUpdateObjectSelection(G,rec->obj);
          break;
        }

        rec->cand_id = TrackerNewCand(I->Tracker,(TrackerRef*)rec);
        TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id,1);
        
        switch(rec->type) {
        case cExecObject:
          TrackerLink(I->Tracker, rec->cand_id, I->all_obj_list_id,1);
          break;
        case cExecSelection:
          TrackerLink(I->Tracker, rec->cand_id, I->all_sel_list_id,1);
          break;
        }
        ListAppend(I->Spec,rec,next,SpecRec);
        ExecutiveAddKey(I,rec);
        ExecutiveInvalidateGroups(G,false);
        ExecutiveInvalidatePanelList(G);
      } else {
        ListElemFree(rec);
      }
    }
    a++;
    if(!ok) {
      incomplete=true;
      ok=true;
    }
  }
  return(!incomplete);
}

static int ExecutiveSetSelections(PyMOLGlobals *G,PyObject *names)
{
  /* must already have objects loaded at this point... */

  int ok=true;
  int a=0,l=0;
  PyObject *cur;
  SpecRec *rec = NULL;
  int extra;
  int incomplete = false;

  if(ok) ok = (names!=NULL);
  if(ok) ok = PyList_Check(names);
  if(ok) l = PyList_Size(names);
  while(ok&&(a<l)) {
    cur = PyList_GetItem(names,a);
    if(cur!=Py_None) { /* skip over None w/o aborting */
      rec=NULL;
      ListElemCalloc(G,rec,SpecRec); 
      rec->next=NULL;

      if(ok) ok = PyList_Check(cur);
      if(ok) ok = PConvPyStrToStr(PyList_GetItem(cur,0),rec->name,sizeof(WordType));
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,1),&rec->type);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,2),&rec->visible);
      if(ok) ok = PConvPyListToIntArrayInPlaceAutoZero(PyList_GetItem(cur,3),
                                                       rec->repOn,cRepCnt);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(cur,4),&extra);
      switch(rec->type) {
      case cExecSelection:
        ok = SelectorFromPyList(G,rec->name,PyList_GetItem(cur,5));
        break;
      }
      ListElemFree(rec);
    }
    a++;
    if(!ok) {
      incomplete=true;
      ok=true;
    }
  }
  return(!incomplete);
}

static PyObject *ExecutiveGetExecSeleAsPyList(PyMOLGlobals *G,SpecRec *rec)
{
  PyObject *result = NULL;
  int sele;

  sele = SelectorIndexByName(G,rec->name);
  if(sele>=0) {
    result = PyList_New(7);
    PyList_SetItem(result,0,PyString_FromString(rec->name));
    PyList_SetItem(result,1,PyInt_FromLong(cExecSelection));
    PyList_SetItem(result,2,PyInt_FromLong(rec->visible));
    PyList_SetItem(result,3,PConvIntArrayToPyList(rec->repOn,cRepCnt));
    PyList_SetItem(result,4,PyInt_FromLong(-1));
    PyList_SetItem(result,5,SelectorAsPyList(G,sele));
    PyList_SetItem(result,6,PyString_FromString(rec->group_name));
    
  }
  return(PConvAutoNone(result));
}

static PyObject *ExecutiveGetNamedEntries(PyMOLGlobals *G,int list_id,int partial)
{
  register CExecutive *I = G->Executive;  
  CTracker *I_Tracker= I->Tracker;
  PyObject *result = NULL;
  int count=0,total_count=0;
  int iter_id = 0;
  SpecRec *rec = NULL, *list_rec = NULL;

  SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);

  if(list_id) {
    total_count = TrackerGetNCandForList(I_Tracker,list_id);
    iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  } else {
    total_count = ExecutiveCountNames(G);
  }
  result = PyList_New(total_count);

  /* critical reliance on short-circuit behavior */

  while( (iter_id && TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&list_rec)) ||
         ((!iter_id) && ListIterate(I->Spec,rec,next))) {
    
    if(list_id)
      rec = list_rec;
    if(count>=total_count)
      break;
    if(rec) {
      switch(rec->type) {
      case cExecObject:
        PyList_SetItem(result,count,
                       ExecutiveGetExecObjectAsPyList(G,rec));
        break;
      case cExecSelection:
        if(!partial) {
          PyList_SetItem(result,count,
                         ExecutiveGetExecSeleAsPyList(G,rec));
        } else { 
          /* cannot currently save selections in partial sessions */
          PyList_SetItem(result,count,PConvAutoNone(NULL));
        }
        break;
      default:
        PyList_SetItem(result,count,PConvAutoNone(NULL));
        break;
      }
    } else {
      PyList_SetItem(result,count,PConvAutoNone(NULL));
    }
    count++;
  }

  while(count<total_count) { /* insure that all members of outgoing list are defined */
    PyList_SetItem(result,count,PConvAutoNone(NULL));
    count++;
  }
  
  if(iter_id) {
    TrackerDelIter(I_Tracker, iter_id);
  }
  return(PConvAutoNone(result));
}


#endif

#ifndef _PYMOL_NOPY
#ifdef PYMOL_EVAL
#include "ExecutiveEvalMessage.h"
#endif
#endif

int ExecutiveGetSession(PyMOLGlobals *G,PyObject *dict,char *names,int partial,int quiet)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok=true;
  int list_id=0;
  SceneViewType sv;
  PyObject *tmp;
  
  if(names && names[0]) {
    list_id = ExecutiveGetNamesListFromPattern(G,names,true,cExecExpandKeepGroups);
  } 

  tmp = PyInt_FromLong(_PyMOL_VERSION_int);
  PyDict_SetItemString(dict,"version",tmp);
  Py_XDECREF(tmp);

  tmp = ExecutiveGetNamedEntries(G,list_id,partial);
  PyDict_SetItemString(dict,"names",tmp);
  Py_XDECREF(tmp);

  tmp = ColorAsPyList(G);
  PyDict_SetItemString(dict,"colors",tmp);
  Py_XDECREF(tmp);
            
  tmp = ColorExtAsPyList(G); 
  PyDict_SetItemString(dict,"color_ext",tmp);
  Py_XDECREF(tmp);

  tmp = SettingUniqueAsPyList(G);
  PyDict_SetItemString(dict,"unique_settings",tmp);
  Py_XDECREF(tmp);

  if(partial) { /* mark this as a partial session */

    PyDict_SetItemString(dict,"partial",PConvAutoNone(Py_None));

  } else { 

    /* none of the following information is saved in partial sessions */

    tmp = SelectorSecretsAsPyList(G);
    PyDict_SetItemString(dict,"selector_secrets",tmp);
    Py_XDECREF(tmp);
    
    tmp = SettingGetGlobalsAsPyList(G);
    PyDict_SetItemString(dict,"settings",tmp);
    Py_XDECREF(tmp);

    SceneGetView(G,sv);
    tmp = PConvFloatArrayToPyList(sv,cSceneViewSize);
    PyDict_SetItemString(dict,"view",tmp);
    Py_XDECREF(tmp);
    
    tmp = MovieAsPyList(G);
    PyDict_SetItemString(dict,"movie",tmp);
    Py_XDECREF(tmp);
    
    tmp = EditorAsPyList(G);
    PyDict_SetItemString(dict,"editor",tmp);
    Py_XDECREF(tmp);
    
#ifndef _PYMOL_NO_MAIN
    tmp = MainAsPyList();
    PyDict_SetItemString(dict,"main",tmp);
    Py_XDECREF(tmp);
#endif

#ifdef PYMOL_EVAL
    ExecutiveEvalMessage(G,dict);
#endif

  }
  
  if(Feedback(G,FB_Executive,FB_Errors)) {
    if(PyErr_Occurred()) {
      PRINTF
        " ExecutiveGetSession: a Python error occured during creation of the session object:\n"
        ENDF(G);
      PyErr_Print();
    }
  }

  return(ok);
#endif

}

#ifndef _PYMOL_NOPY
static void ExecutiveMigrateSession(PyMOLGlobals *G,int session_version)
{
  if(session_version<100) {
    /* migrate lighting model */
    SettingSetGlobal_f(G,cSetting_direct,
                       1.8*SettingGetGlobal_f(G,cSetting_direct));
    SettingSetGlobal_f(G,cSetting_reflect,
                       0.5*SettingGetGlobal_f(G,cSetting_reflect));
    SettingSetGlobal_f(G,cSetting_ambient,
                       1.166*SettingGetGlobal_f(G,cSetting_ambient));
    SettingSetGlobal_f(G,cSetting_gamma,
                       0.769*SettingGetGlobal_f(G,cSetting_gamma));

    /* try best to meet existing expectations with existing sessions */
    SettingSetGlobal_f(G,cSetting_ray_legacy_lighting, 1.0F);

    /* force use of movie_delay in preference to movie_fps */

    SettingSetGlobal_f(G,cSetting_movie_fps, 0.0F);

    /* and labels */
    SettingSetGlobal_i(G,cSetting_label_digits, 2);
    SettingSetGlobal_3f(G,cSetting_label_position, 1.0F,1.0F,0.0F);
  }
  if(session_version<99) {
    SettingSetGlobal_f(G,cSetting_cartoon_ladder_mode,0);
    SettingSetGlobal_f(G,cSetting_cartoon_tube_cap,0);
    SettingSetGlobal_f(G,cSetting_cartoon_nucleic_acid_mode,1);
    {
      float old_sulfur[3] = {1.0, 0.5, 0.0};
      ColorDef(G,"sulfur",old_sulfur,0,true);  
    }
  }
  if(session_version<98) { /* produce expected rendering quality & performance with old sessions */
    SettingSetGlobal_b(G,cSetting_ray_orthoscopic,1);
  }
  if(session_version<96) {
    SettingSetGlobal_f(G,cSetting_ray_transparency_contrast, 1.0F);
  }
  if(session_version<95) {

    { /* adjust fog to reflect current importance of seeing to the Z-slab center w/o fog */
      
      float fog_start = SettingGetGlobal_f(G,cSetting_fog_start);
      float ray_trace_fog_start = SettingGetGlobal_f(G,cSetting_ray_trace_fog_start);
      if((fog_start==0.40F)||(fog_start==0.35F)||(fog_start==0.30F)) {
        SettingSetGlobal_f(G,cSetting_fog_start,0.45F);
      }
      if((ray_trace_fog_start==0.45F)||(ray_trace_fog_start==0.40F)||(ray_trace_fog_start==0.35F)) {
        SettingSetGlobal_f(G,cSetting_ray_trace_fog_start,0.50F);
      }

    }

    { /* adjust GUI width */

      int gui_width = SettingGetGlobal_i(G,cSetting_internal_gui_width);

      if(gui_width==160) {
        SettingSetGlobal_i(G,cSetting_internal_gui_width,220);
      }
    }

    { /* enable antialiasing */

      int antialias = SettingGetGlobal_i(G,cSetting_antialias);

      if(antialias==0) {
        SettingSetGlobal_i(G,cSetting_antialias,1);
      }
      
    }
  }
}
#endif

int ExecutiveSetSession(PyMOLGlobals *G,PyObject *session,
                        int partial_restore,int quiet)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok=true;
  int incomplete = false;
  PyObject *tmp;
  SceneViewType sv;
  int version=-1;
  int migrate_sessions = SettingGetGlobal_b(G,cSetting_session_migration);
  char active[WordLength] = "";
  int  have_active = false;
  int partial_session = false;

  if(!partial_restore) { /* if user has requested partial restore */
    ExecutiveDelete(G,"all");
    ColorReset(G);
  }

  if(ok) ok = PyDict_Check(session);

  if(ok) { /* if session is partial, then don't error about missing stuff */
    tmp = PyDict_GetItemString(session,"partial"); 
    if(tmp) {
      partial_session = true;
    }
  }

  if(ok) {
    tmp = PyDict_GetItemString(session,"version");
    if(tmp) {
      ok = PConvPyIntToInt(tmp,&version);
      if(ok) {
        if(version>_PyMOL_VERSION_int) {
          if(!quiet) {
            PRINTFB(G,FB_Executive,FB_Errors)
              "Warning: This session was created with a newer version of PyMOL (%1.2f).\n",
              version/100.
              ENDFB(G);
            if(SettingGet(G,cSetting_session_version_check)) {
              PRINTFB(G,FB_Executive,FB_Errors)
                "Error: Please update first -- see http://www.pymol.org\n"
                ENDFB(G);
              ok=false;
            } else {
              PRINTFB(G,FB_Executive,FB_Errors)
                "Warning: Some content may not load completely.\n"
                ENDFB(G);
            }
          }
        } else {
          if(!quiet) {
            PRINTFB(G,FB_Executive,FB_Details)          
              " Executive: Loading version %1.2f session...\n",
              version/100.0
              ENDFB(G);
          }
        }
      }
    }
  }

#ifndef PYMOL_EVAL
  if(ok) {
    tmp = PyDict_GetItemString(session,"eval_nag");
    if(tmp) {
      ok = PyString_Check(tmp);
      if(ok) {
        char *st = PyString_AsString(tmp);
        if(st) {
          if(Feedback(G,FB_Nag,FB_Warnings)) {
            OrthoAddOutput(G,st);
          }
        }
      }
    }
  }
#endif

  if(ok) { /* colors must be restored before settings and objects */
    tmp = PyDict_GetItemString(session,"colors");
    if(tmp) {
      ok = ColorFromPyList(G,tmp,partial_restore);
    } 

    if(tmp||(!partial_restore)) { /* ignore missing if partial restore */
    
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after colors.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }

  if(ok) {
    tmp = PyDict_GetItemString(session,"color_ext");
    if(tmp) {
      ok = ColorExtFromPyList(G,tmp,partial_restore);
    }

    if(tmp||(!partial_session)) { /* ignore missing if partial restore */    
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after color_ext.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"settings");
    if(tmp) {
      ok = SettingSetGlobalsFromPyList(G,tmp);
    }

    if(tmp||(!(partial_restore|partial_session))) { /* ignore missing if partial restore */        
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after settings.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"unique_settings");
    if(tmp) {
      ok = SettingUniqueFromPyList(G,tmp,partial_restore);
    }
    
    if(tmp||(!partial_session)) { /* ignore missing if partial restore */    
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after settings.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }
  if(ok) {
    tmp = PyDict_GetItemString(session,"names");
    if(tmp) {
      if(ok) ok=ExecutiveSetNamedEntries(G,tmp,version,partial_restore,partial_session);
      if(!(partial_restore||partial_session)) {
        if(ok) ok=ExecutiveSetSelections(G,tmp);
        if(ok) have_active = ExecutiveGetActiveSeleName(G,active,false,false);
      }
    }
    if(PyErr_Occurred()) {
      PyErr_Print();  
      ok=false;
    }
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "ExectiveSetSession-Error: after names.\n"
        ENDFB(G);
    }
    if(!ok) {
      incomplete = true;
      ok=true; /* keep trying...don't give up */
    }
  }
  if(ok&&!(partial_restore)) {
    tmp = PyDict_GetItemString(session,"selector_secrets");
    if(tmp) {
      if(ok) ok=SelectorSecretsFromPyList(G,tmp);
    }
  
    if(tmp||(!(partial_restore|partial_session))) { /* ignore missing if partial restore */      
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after selector secrets.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }  
  if(ok&&!(partial_restore)) {
    tmp = PyDict_GetItemString(session,"view");
    if(tmp) {
      ok = PConvPyListToFloatArrayInPlace(tmp,sv,cSceneViewSize);
    }
    if(tmp||(!(partial_restore|partial_session))) { /* ignore missing if partial restore */    
      if(ok) SceneSetView(G,sv,true,0,0);
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after view.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }
  if(ok&&!(partial_restore)) {
    int warning;
    tmp = PyDict_GetItemString(session,"movie");
    if(tmp) {
      ok = MovieFromPyList(G,tmp,&warning);
    }
    if(tmp||(!(partial_restore|partial_session))) { /* ignore missing if partial restore */    
      
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after movie.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }
  
  if(ok&&!(partial_restore)) {
    tmp = PyDict_GetItemString(session,"editor");
    if(tmp) {
      ok = EditorFromPyList(G,tmp);
    }
    if(tmp||(!(partial_restore|partial_session))) { /* ignore missing if partial restore */    
      
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after editor.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }
  if(ok) { /* update mouse in GUI */
    PParse(G,"cmd.mouse(quiet=1)");
    if(!G->Option->presentation)
      PParse(G,"viewport"); /* refresh window/internal_gui status */
  }
 
#ifndef _PYMOL_NO_MAIN
  if(ok) {
    tmp = PyDict_GetItemString(session,"main");
    if(tmp) {
      ok = MainFromPyList(tmp);  /* main just stores the viewport size */
    }
    if(tmp||(!(partial_restore|partial_session))) { /* ignore missing if partial restore */    
      if(PyErr_Occurred()) {
        PyErr_Print();  
        ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "ExectiveSetSession-Error: after main.\n"
          ENDFB(G);
      }
      if(!ok) {
        incomplete = true;
        ok=true; /* keep trying...don't give up */
      }
    }
  }
#endif

  if(ok&&migrate_sessions) { /* migrate sessions */
    tmp = PyDict_GetItemString(session,"version");
    if(tmp) {
      ok = PConvPyIntToInt(tmp,&version);
      if(ok) {
        ExecutiveMigrateSession(G,version);
      }
    }
  }
  if(ok) {
    if(have_active)
      ExecutiveSetObjVisib(G,active,true,false);      
  }
  if(incomplete) {
    PRINTFB(G,FB_Executive,FB_Warnings)
      "ExectiveSetSession-Warning: restore may be incomplete.\n"
      ENDFB(G);
  }
  return(ok);
#endif
}

#define ExecScrollBarMargin 1
#define ExecScrollBarWidth 13

void ExecutiveObjMolSeleOp(PyMOLGlobals *G,int sele,ObjectMoleculeOpRec *op);

static CObject **ExecutiveSeleToObjectVLA(PyMOLGlobals *G,char *s1)
{
  /* return VLA containing list of atoms references by selection */

  CObject **result = NULL;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  CObject *obj=NULL;
  int n = 0;
  ObjectMoleculeOpRec op2;
  int sele;

  result = VLAlloc(CObject*,50);
  if(WordMatch(G,s1,cKeywordAll,true)) {
    /* all objects */
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject)
          {
            VLACheck(result,CObject*,n);
            result[n]=rec->obj;
            n++;
          }
      }
  } else {
    sele = SelectorIndexByName(G,s1);
    if(sele>0) {
      ObjectMoleculeOpRecInit(&op2);
      op2.code=OMOP_GetObjects;
      op2.obj1VLA=(ObjectMolecule**)result;
      op2.i1=0;
      ExecutiveObjMolSeleOp(G,sele,&op2);
      n = op2.i1;
      result = (CObject**)op2.obj1VLA;
    } else {
      obj = ExecutiveFindObjectByName(G,s1);
      if(obj) {
        VLACheck(result,CObject*,n);
        result[n]=obj;
        n++;
      }
    }
  }
  VLASize(result,CObject*,n);
  return(result);
}

int ExecutiveGetCrystal(PyMOLGlobals *G,char *sele,float *a,float *b,float *c,
                        float *alpha,float *beta,float *gamma,
                        char *sgroup,int *defined)
{
  int ok=true;

  ObjectMolecule *objMol;
  int sele0;
  sele0 = SelectorIndexByName(G,sele);
  *defined = false;
  if(sele0<0) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Error: invalid selection.\n"
      ENDFB(G);
    ok=false;
  } else {
    objMol = SelectorGetSingleObjectMolecule(G,sele0);
    if(!objMol) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "Error: selection must refer to exactly one object.\n"
        ENDFB(G);
      ok=false;
    } else {
      if(objMol->Symmetry&&objMol->Symmetry->Crystal) {
        *a = objMol->Symmetry->Crystal->Dim[0];
        *b = objMol->Symmetry->Crystal->Dim[1];
        *c = objMol->Symmetry->Crystal->Dim[2];
        *alpha = objMol->Symmetry->Crystal->Angle[0];
        *beta = objMol->Symmetry->Crystal->Angle[1];
        *gamma = objMol->Symmetry->Crystal->Angle[2];
        UtilNCopy(sgroup, objMol->Symmetry->SpaceGroup,sizeof(WordType));
        *defined = true;
      }
    }
  }
  return(ok);
}

int ExecutiveSetCrystal(PyMOLGlobals *G,char *sele,float a,float b,float c,
                        float alpha,float beta,float gamma,char *sgroup)
{
  CObject **objVLA = NULL;
  CObject *obj;
  ObjectMolecule *objMol;
  /*  ObjectMap *objMap;*/
  int ok=true;
  CSymmetry *symmetry = NULL;
  CCrystal *crystal = NULL;
  int n_obj;
  int i;

  objVLA = ExecutiveSeleToObjectVLA(G,sele);
  n_obj = VLAGetSize(objVLA);
  if(n_obj) {
    for(i=0;i<n_obj;i++) {
      obj = objVLA[i];
      switch(obj->type) {
      case cObjectMolecule:
        if(!symmetry) {
          symmetry=SymmetryNew(G);          
          symmetry->Crystal->Dim[0]=a;
          symmetry->Crystal->Dim[1]=b;
          symmetry->Crystal->Dim[2]=c;
          symmetry->Crystal->Angle[0]=alpha;
          symmetry->Crystal->Angle[1]=beta;
          symmetry->Crystal->Angle[2]=gamma;
          UtilNCopy(symmetry->SpaceGroup,sgroup,sizeof(WordType));
          SymmetryAttemptGeneration(symmetry,false);
        }
        objMol = (ObjectMolecule*)obj;
        if(symmetry) {
          if(objMol->Symmetry)
            SymmetryFree(objMol->Symmetry);
          objMol->Symmetry = SymmetryCopy(symmetry);
        }
        break;
        /* not currently supported 
           case cObjectMap:
        
           if(!crystal) {
           crystal = CrystalNew(G);
           crystal->Dim[0]=a;
           crystal->Dim[1]=b;
           crystal->Dim[2]=c;
           crystal->Angle[0]=alpha;
           crystal->Angle[1]=beta;
           crystal->Angle[2]=gamma;
           CrystalUpdate(crystal);
           }
           if(crystal) {
           objMap = (ObjectMap*)obj;
           if(objMap->Crystal) {
           CrystalFree(objMap->Crystal);
           }
           objMap->Crystal=CrystalCopy(crystal);
           }
           break;
        */

      }
    }
  } else {
    ok=false;
    PRINTFB(G,FB_Executive,FB_Errors)
      " ExecutiveSetCrystal: no object selected\n"
      ENDFB(G);
  }
  if(crystal)
    CrystalFree(crystal);
  if(symmetry)
    SymmetryFree(symmetry);
  VLAFreeP(objVLA);
  return(ok);
}

int ExecutiveSmooth(PyMOLGlobals *G,char *name,int cycles,
                    int window,int first, int last, int ends,
                    int quiet)
{

  int sele = -1;
  ObjectMoleculeOpRec op;
  int state;
  int n_state;
  float *coord0=NULL,*coord1=NULL;
  int *flag0=NULL,*flag1=NULL;
  int a,b,c,d,st,cnt;
  float i_cnt;
  int n_atom;
  int ok=true;
  int backward;
  int forward;
  int range,offset;
  int end_skip=0;
  float *v0,*v1;
  float sum[3];
  int loop = false;
  /*  WordType all = "_all";*/

  PRINTFD(G,FB_Executive)
    " ExecutiveSmooth: entered %s,%d,%d,%d,%d,%d\n",name,cycles,first,last,window,ends
    ENDFD;

  sele=SelectorIndexByName(G,name);

  if(sele>=0) {
    int max_state = ExecutiveCountStates(G,name)-1;
    if(last<0) 
      last = max_state;
    if(first<0)
      first = 0;
    if(last<first) {
      state=last;
      last=first;
      first=state;
    }
    if(last>max_state)
      last = max_state;

    n_state=last-first+1;

    backward=window/2;
    forward=window/2;

    if((forward-backward)==(window+1))
      forward--; /* even sizes window */
    
    switch(ends) {
    case 0:
      end_skip = 1;
      break;
    case 1:
      end_skip = 0;
      break;
    case 2:
      end_skip = backward;
      break;
    case 3: /* cyclic averaging */
      end_skip = 0;
      loop = true;
      break;
    default:
      end_skip = 0;
      break;
    }

    if(ends) {
      range = (last-first)+1;
      offset = 0;
    } else {
      range = (last-end_skip)-(first+end_skip)+1;
      offset = end_skip;
    }
    
    PRINTFD(G,FB_Executive)
      " ExecutiveSmooth: first %d last %d n_state %d backward %d forward %d range %d\n",
      first,last,n_state,backward,forward,range
      ENDFD;

    if(n_state>=window) {
      
      /* determine storage req */
      ObjectMoleculeOpRecInit(&op);
      op.code = OMOP_CountAtoms;
      op.i1=0;
      ExecutiveObjMolSeleOp(G,sele,&op);
      n_atom=op.i1;
      if(n_atom) {
        /* allocate storage */
        coord0 = Alloc(float,3*n_atom*n_state);
        coord1 = Alloc(float,3*n_atom*n_state);
        flag0 = Alloc(int,n_atom*n_state);
        flag1 = Alloc(int,n_atom*n_state);
        
        /* clear the arrays */
        
        UtilZeroMem(coord0,sizeof(float)*3*n_atom*n_state);
        UtilZeroMem(flag0,sizeof(int)*n_atom*n_state);
        
        /* get the data */
        if(!quiet) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Smooth: copying coordinates to temporary arrays..\n"
            ENDFB(G);
        }
        op.code = OMOP_CSetIdxGetAndFlag;
        op.i1 = n_atom; 
        op.i2 = 0;
        op.cs1 = first;
        op.cs2 = last;
        op.vv1 = coord0;
        op.ii1 = flag0;
        op.nvv1 = 0;          
        ExecutiveObjMolSeleOp(G,sele,&op);    
        
        PRINTFD(G,FB_Executive)  
          " ExecutiveSmooth: got %d %d\n",op.i2,op.nvv1
          ENDFD;
        
        UtilZeroMem(coord1,sizeof(float)*3*n_atom*n_state);
        UtilZeroMem(flag1,sizeof(int)*n_atom*n_state);
        
        for(a=0;a<cycles;a++) {                
          if(!quiet) {
            PRINTFB(G,FB_Executive,FB_Actions)
              " Smooth: smoothing (pass %d)...\n",a+1
              ENDFB(G);
          }
          for(b=0;b<range;b++) {
            for(c=0;c<n_atom;c++) {
              zero3f(sum);
              cnt = 0;
              for(d=-backward;d<=forward;d++) {
                st = b + offset + d;
                if(loop) {
                  if(st<0) {
                    st=n_state+st;
                  } else if(st>=n_state) {
                    st=st-n_state;
                  }
                } else {
                  if(st<0) {
                    st=0;
                  } else if(st>=n_state) {
                    st=n_state-1;
                  }
                }

                /*if(c==0) printf("averaging from slot %d\n",st);*/
                cnt+=flag0[(n_atom*st)+c];
                v0 = coord0 + 3*(n_atom*st+c);
                add3f(sum,v0,sum);
              }
              if(cnt) {
                st = b + offset;
                if((st>=end_skip)&&(st<(n_state-end_skip))) {
                  /* if(c==0) printf("dumping into slot %d\n",st);*/
                  flag1[(n_atom*st)+c] = flag0[(n_atom*st)+c]; /* don't flag states that weren't originally flagged */
                  i_cnt = 1.0F/cnt;
                  v1 = coord1 + 3*((n_atom*st)+c);
                  scale3f(sum,i_cnt,v1);
                }
              }
            }
          }
          for(b=0;b<range;b++) {
            for(c=0;c<n_atom;c++) {
              st = b + offset;
              if(flag1[(n_atom*st)+c]) {
                v0 = coord0 + 3*((n_atom*st)+c);
                v1 = coord1 + 3*((n_atom*st)+c);
                copy3f(v1,v0);
              }
            }
          }
        }

        if(!quiet) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Smooth: updating coordinates...\n"
            ENDFB(G);
        }
        
        /* set the new coordinates */
        
        op.code = OMOP_CSetIdxSetFlagged;
        op.i1 = n_atom; 
        op.i2 = 0;
        if(ends) {
          op.cs1 = first;
          op.cs2 = last;
          op.vv1 = coord1;
          op.ii1 = flag1;
        } else {
          op.cs1 = first+end_skip;
          op.cs2 = last-end_skip;
          op.vv1 = coord1+(end_skip*3*n_atom);
          op.ii1 = flag1+(end_skip*n_atom);
        }
        op.nvv1 = 0;
        
        ExecutiveObjMolSeleOp(G,sele,&op);      
        PRINTFD(G,FB_Executive)  
          " ExecutiveSmooth: put %d %d\n",op.i2,op.nvv1
          ENDFD;
        
        FreeP(coord0);
        FreeP(coord1);
        FreeP(flag0);
        FreeP(flag1);
      }
    }
  } else {
    PRINTFB(G,FB_Executive,FB_Errors)  
      " ExecutiveSmooth: selection not found\n"
      ENDFB(G);
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveDebug(PyMOLGlobals *G,char *name)
{
  ObjectMolecule *obj;
  ObjectMoleculeBPRec bp;
  int a;

  obj=(ObjectMolecule*)ExecutiveFindObjectByName(G,name);
  if(obj) {
    ObjectMoleculeInitBondPath(obj,&bp);
    ObjectMoleculeGetBondPaths(obj,0,10,&bp);
    for(a=0;a<bp.n_atom;a++) {
      printf("%d %d %d\n",a,bp.list[a],bp.dist[bp.list[a]]);
    }
    
    ObjectMoleculePurgeBondPath(obj,&bp);
  }
  return(1);
}
/*========================================================================*/
int ***ExecutiveGetBondPrint(PyMOLGlobals *G,char *name,int max_bond,int max_type,int *dim)
{
  int ***result = NULL;
  CObject *obj;
  ObjectMolecule *objMol;

  obj=ExecutiveFindObjectByName(G,name);
  if(obj->type==cObjectMolecule) {
    objMol = (ObjectMolecule*)obj;
    result = ObjectMoleculeGetBondPrint(objMol,max_bond,max_type,dim);
  }
  return(result);
}
/*========================================================================*/
#define cMapOperatorMinimum 0
#define cMapOperatorMaximum 1
#define cMapOperatorSum     2
#define cMapOperatorAverage 3
#define cMapOperatorDifference 4
#define cMapOperatorCopy     5
#define cMapOperatorUnique   6

int ExecutiveMapSet(PyMOLGlobals *G,char *name,int operator,char *operands,
                    int target_state,int source_state,int zoom, int quiet)
{
  CExecutive *I = G->Executive;
  CTracker *I_Tracker= I->Tracker;
  int ok=true;
  int isNew = false;
  ObjectMap *target = ExecutiveFindObjectMapByName(G,name);
  ObjectMap *first_operand = NULL;
  int src_state_start=0,src_state_stop=0;
  int list_id = ExecutiveGetNamesListFromPattern(G,operands,true,true);
  
  if(target_state<0) /* if we're targeting all states, 0 is the offset */
    target_state = 0;

  /* first, figure out what the range of input states is */
  
  if(source_state<0) { /* all source states */
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int max_n_state = 0;
    SpecRec *rec;    
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        switch(rec->type) {
        case cExecObject:
          if(rec->obj->type==cObjectMap) {
            ObjectMap *obj =(ObjectMap*)rec->obj;
            if(obj->NState>max_n_state)
              max_n_state = obj->NState; /* count states */
          }
        }
      }
    }
    TrackerDelIter(I_Tracker, iter_id);
    src_state_start = 0;
    src_state_stop = max_n_state;
  } else {
    src_state_start = source_state;
    src_state_stop = source_state+1;
  }

  {
    /* next, find the first operand */
    
    OrthoLineType first_op_st;
    ParseWordCopy(first_op_st,operands,sizeof(OrthoLineType)-1); /* copy the first operand */
    {
      int sub_list_id = ExecutiveGetNamesListFromPattern(G,first_op_st,true,true);
      int sub_iter_id = TrackerNewIter(I_Tracker, 0, sub_list_id);
      SpecRec *rec;
      
      while( TrackerIterNextCandInList(I_Tracker, sub_iter_id, (TrackerRef**)&rec) ) {
        if(rec) {
          switch(rec->type) {
          case cExecObject:
            if(rec->obj->type==cObjectMap) {
              ObjectMap *obj =(ObjectMap*)rec->obj;
              first_operand = obj;
            }
            break;
          }
        }
        if(first_operand) break;
      }
      TrackerDelList(I_Tracker, sub_list_id);
      TrackerDelIter(I_Tracker, sub_iter_id);
    }
  }

  {

    /* okay, next thing we need to worry about is where we're putting the data.

    Case 1. If the map already exists, then we'll use the existing map points for storing
    the result.

    Case 2. If the operation implies a copy of existing map geometry, then we'll create that
    copy first before performing the calulation.

    Case 3. If the operation implies a new map geometry, then we need to compute that geometry and
    create the map.

    */

    if(!target) { /* target map doesn't exist...*/
      int need_union_geometry = false;
      int need_first_geometry = false;
      switch(operator) {
      case cMapOperatorSum:
      case cMapOperatorAverage:
      case cMapOperatorMinimum:
      case cMapOperatorMaximum:
      case cMapOperatorDifference:
        need_union_geometry = true;
        break;
      case cMapOperatorUnique:
      case cMapOperatorCopy:
        need_first_geometry = true;
        break;
      }
      
      if(need_union_geometry) {
        int src_state,trg_state;
        ObjectMapDesc desc;
        target = ObjectMapNew(G);
        
        ObjectSetName((CObject*)target,name);
        isNew = true;

        for(src_state=src_state_start;src_state<src_state_stop;src_state++) {
          trg_state = src_state + target_state;
          desc.mode = cObjectMap_OrthoMinMaxGrid; /* Orthorhombic: min, max, 
                                                     spacing, 
                                                     centered over range  */
          desc.init_mode = 0; /* zeros */
          desc.Grid[0] = 1.0F;
          desc.Grid[1] = 1.0F;
          desc.Grid[2] = 1.0F;

          {
            int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
            float grid_sum[3] = {0.0F,0.0F,0.0F};
            int grid_count = 0;
            int first_extent = true;

            SpecRec *rec;    
            while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
              if(rec) {
                /* compute an average grid and get the effective "union" extent */
                switch(rec->type) {
                case cExecObject:
                  if(rec->obj->type==cObjectMap) {
                    
                    ObjectMap *obj =(ObjectMap*)rec->obj;
                    ObjectMapState *ms = obj->State + src_state;
                    if(src_state<obj->NState) {
                      if(ms->Active) {
                        if(first_extent) {
                          copy3f(ms->ExtentMin,desc.MinCorner);
                          copy3f(ms->ExtentMax,desc.MaxCorner);
                          first_extent = false;
                        } else {
                          int b;
                          for(b=0;b<3;b++) {
                            if(ms->ExtentMin[b]<desc.MinCorner[b])
                              desc.MinCorner[b]=ms->ExtentMin[b];
                            if(ms->ExtentMax[b]>desc.MaxCorner[b])
                              desc.MaxCorner[b]=ms->ExtentMax[b];
                          }
                        }
                        switch(ms->MapSource) {
                        case cMapSourcePHI:
                        case cMapSourceFLD:
                        case cMapSourceDesc:
                        case cMapSourceChempyBrick:
                        case cMapSourceVMDPlugin:
                          {
                            int b;
                            for(b=0;b<3;b++) {
                              grid_sum[b] += ms->Grid[b];
                            }
                          }
                          grid_count++;
                          break;
                          /* other map state types not currently handled... */
                        }
                      }
                    }
                  }
                }
              }
            }
            TrackerDelIter(I_Tracker, iter_id);
            if(grid_count) {
              int b;
              for(b=0;b<3;b++) {
                desc.Grid[b] = grid_sum[b]/grid_count;
              }
            }
            if(!first_extent) {
              float tmp[3];
              scale3f(desc.Grid,0.5,tmp);
              add3f(desc.Grid,desc.MaxCorner,desc.MaxCorner);
              subtract3f(desc.MinCorner,desc.Grid,desc.MinCorner);
              ObjectMapNewStateFromDesc(G,target,&desc,trg_state,quiet);
              if(trg_state>=target->NState)
                target->NState = trg_state+1;
              target->State[trg_state].Active=true;
            }
          }
        }
        /* need union geometry */
      } else if(need_first_geometry) {
        if(first_operand) {
          if( ObjectMapNewCopy(G,first_operand, &target, source_state, target_state )) {
            if(target) {
              ObjectSetName((CObject*)target,name);
              isNew = true;
            }
          }
        }
      }
    }
  }

  if(!target) {
    ok=false;
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: cannot find or construct target map.\n"   
      ENDFB(G);
  } 

  /* now do the actual operation */
  
  if(ok&&target) {
    int src_state;
    for(src_state=src_state_start;src_state<src_state_stop;src_state++) {
      int trg_state = src_state + target_state;
      ObjectMapState *ms;
      VLACheck(target->State,ObjectMapState,trg_state);
      
      ms = target->State + target_state;    
      if(ms->Active) {
        int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
        int n_pnt = (ms->Field->points->size / ms->Field->points->base_size)/3;
        float *pnt = (float*)ms->Field->points->data;
        float *r_value = Alloc(float, n_pnt);
        float *l_value = Calloc(float, n_pnt);
        int *present = Calloc(int, n_pnt);
        int *inside = Alloc(int, n_pnt);
        SpecRec *rec;
        
        while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
          if(rec) {
            if(rec->type == cExecObject) {
              if(rec->obj->type==cObjectMap) {
                ObjectMap *obj =(ObjectMap*)rec->obj;
                ObjectMapInterpolate(obj, src_state, pnt, r_value, inside, n_pnt);
                {
                  register int a;
                  register float *rv = r_value;
                  register float *lv = l_value;
                  register int *flg = inside;
                  register int *pre = present;
                
                  switch(operator) {
                  case cMapOperatorCopy:
                    for(a=0;a<n_pnt;a++) {
                      if(flg) {
                        *lv = *rv;
                      }
                      rv++; lv++; flg++;
                    }
                    break;
                  case cMapOperatorMinimum:
                    for(a=0;a<n_pnt;a++) {
                      if(flg) {
                        if(*pre) { 
                          if(*lv>*rv)
                            *lv = *rv;
                        } else { /* first map */
                          *pre = 1;
                          *lv = *rv;
                        }
                      
                      }
                      rv++; lv++; flg++; pre++;
                    }
                    break;
                  case cMapOperatorMaximum:
                    for(a=0;a<n_pnt;a++) {
                      if(flg) {
                        if(*pre) {
                          if(*lv<*rv)
                            *lv = *rv;
                        } else { /* first map */
                          *pre = 1;
                          *lv = *rv;
                        }
                      }
                      rv++; lv++; flg++; pre++;
                    }
                    break;
                  case cMapOperatorSum:
                    for(a=0;a<n_pnt;a++) {
                      if(flg) {
                        *lv += *rv;
                      }
                      rv++; lv++; flg++;
                    }
                    break;
                  case cMapOperatorAverage:
                    for(a=0;a<n_pnt;a++) {
                      if(flg) {
                        *lv += *rv;
                      }
                      (*pre)++;
                      rv++; lv++; flg++; pre++;
                    }
                    break;
                  case cMapOperatorDifference:
                    if(obj!=first_operand) {
                      for(a=0;a<n_pnt;a++) {
                        if(flg) {
                          *lv -= *rv;
                        }
                        rv++; lv++; flg++;
                      }
                    } else {
                      for(a=0;a<n_pnt;a++) {
                        if(flg) {
                          *lv += *rv;
                        }
                        rv++; lv++; flg++;
                      }
                    }
                    break;
                  case cMapOperatorUnique:
                    if(obj!=first_operand) {
                      for(a=0;a<n_pnt;a++) {
                        if(flg) {
                          *lv -= *rv;
                        }
                        rv++; lv++; flg++;
                      }
                    } else {
                      for(a=0;a<n_pnt;a++) {
                        if(flg) {
                          *lv += *rv;
                        }
                        rv++; lv++; flg++;
                      }
                    }
                    
                    break;
                  }
                }
              }
            }
          }
        }

                  
        
        {
          register int a;
          register float *lv = l_value;
          register int *pre = present;
          
          switch(operator) {
          case cMapOperatorUnique:
            lv = l_value;
            for(a=0;a<n_pnt;a++) {
              if(*lv<0.0F)
                *lv = 0.0F;
              lv++;
            }
            break;
          case cMapOperatorAverage:
            lv = l_value;
            pre = present;
            for(a=0;a<n_pnt;a++) {
              if(*pre)
                *lv /= *pre;
              lv++; pre++;
            }
          }
        }

        /* copy after calculation so that operand can include target */

        memcpy(ms->Field->data->data,l_value,n_pnt*sizeof(float));
     
        FreeP(present);
        FreeP(l_value);
        FreeP(r_value);
        FreeP(inside);
        TrackerDelIter(I_Tracker, iter_id);
      }
    }
  }


  /* and finally, update */

  if(target) {
    ObjectMapUpdateExtents(target);
    if(isNew) {
      ExecutiveManageObject(G,&target->Obj, -1, quiet);
    } else {
      ExecutiveDoZoom(G,&target->Obj,false,zoom,true);
    }
    SceneChanged(G);
  }
  TrackerDelList(I_Tracker, list_id);

  return ok;
}

int ExecutiveMapNew(PyMOLGlobals *G,char *name,int type,float *grid,
                    char *sele,float buffer,
                    float *minCorner,
                    float *maxCorner,int state,int have_corners,
                    int quiet,int zoom,int normalize,float clamp_floor, float clamp_ceiling)
{
  CObject *origObj=NULL;
  ObjectMap *objMap;
  ObjectMapState *ms = NULL;
  int a;
  float v[3];
  ObjectMapDesc _md,*md;
  int ok = true;
  int sele0 = SelectorIndexByName(G,sele);
  int isNew=true;
  int n_state;
  int valid_extent=false;
  int st;
  int st_once_flag=true;
  int n_st;
  int extent_state;
  int clamp_flag = (clamp_floor <= clamp_ceiling);
  
  md=&_md;

  if((state==-2)||(state==-3)) /* TO DO: support per-object states */
    state=SceneGetState(G);

  /* remove object if it already exists */

  origObj=ExecutiveFindObjectByName(G,name);

  if(origObj) {
    if(origObj->type!=cObjectMap) {
      ExecutiveDelete(G,origObj->Name);
    } else {
      isNew=false;
    }
  }

  n_st = ExecutiveCountStates(G,NULL);

  for(st=0;st<n_st;st++) {
    if(state==-1) st_once_flag=false; /* each state, separate map, separate extent */
    if(!st_once_flag) state=st;
    extent_state = state;
    if(state<=-3) extent_state = -1;
    if(strlen(sele)&&(!have_corners)) {
      valid_extent = ExecutiveGetExtent(G,sele,md->MinCorner,
                                        md->MaxCorner,true,extent_state,false);
      /* TODO restrict to state */
    } else {
      valid_extent = 1;
      copy3f(minCorner,md->MinCorner);
      copy3f(maxCorner,md->MaxCorner);
    }
    copy3f(grid,md->Grid);

    subtract3f(md->MaxCorner,md->MinCorner,v);
    for(a=0;a<3;a++) { if(v[a]<0.0) swap1f(md->MaxCorner+a,md->MinCorner+a); };
    subtract3f(md->MaxCorner,md->MinCorner,v);

    if(buffer!=0.0F) {
      for(a=0;a<3;a++) {
        md->MinCorner[a]-=buffer;
        md->MaxCorner[a]+=buffer;
      }
    }
    md->mode = cObjectMap_OrthoMinMaxGrid;
    md->init_mode=-1; /* no initialization */

    /* validate grid */
    for(a=0;a<3;a++) 
      if(md->Grid[a]<=R_SMALL8) md->Grid[a]=R_SMALL8;

    if(ok) {
      if(isNew)
        objMap = ObjectMapNew(G);
      else
        objMap = (ObjectMap*)origObj;
      if(objMap) {
        int once_flag=true;
        n_state = SelectorCountStates(G,sele0);
        if(valid_extent)
          for(a=0;a<n_state;a++) {
            if(state==-5) once_flag=false; /* api: state=-4 = each state, separate map, shared extent */
            if(state==-4) state=-1; /* api: state=-3 all states, but one map */
            if(!once_flag) state=a;
            ms = ObjectMapNewStateFromDesc(G,objMap,md,state,quiet);
            if(!ms)
              ok=false;
          
            if(ok&&ms) {
            
              switch(type) {
              case 0: /* vdw */
                SelectorMapMaskVDW(G,sele0,ms,0.0F,state);
                break;
              case 1: /* coulomb */
                SelectorMapCoulomb(G,sele0,ms,0.0F,state,false,
                                   false,1.0F);
                break;
              case 2: /* gaussian */
                SelectorMapGaussian(G,sele0,ms,0.0F,state,normalize,false,quiet);
                break;
              case 3: /* coulomb_neutral */
                SelectorMapCoulomb(G,sele0,ms,0.0F,state,true, false,1.0F);
                break;
              case 4: /* coulomb_local */
                SelectorMapCoulomb(G,sele0,ms,
                                   SettingGetGlobal_f(G,cSetting_coulomb_cutoff),state,false,
                                   true, 2.0F);
                break;
              case 5: /* gaussian_max */
                SelectorMapGaussian(G,sele0,ms,0.0F,state,normalize,true,quiet);
                break;
              }
              if(!ms->Active)
                ObjectMapStatePurge(G,ms);
              else if(clamp_flag) {
                ObjectMapStateClamp(ms, clamp_floor, clamp_ceiling);
              }
            }
            if(once_flag) break;
          }

        ObjectSetName((CObject*)objMap,name);
        ObjectMapUpdateExtents(objMap);
        if(isNew) {
          ExecutiveManageObject(G,(CObject*)objMap,-1,quiet);
        } else {
          ExecutiveDoZoom(G,(CObject*)objMap,false,zoom,true);
        }
        isNew=false;
        origObj = (CObject*)objMap;
      }
      SceneChanged(G);
    }
    if(st_once_flag)
      break;
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSculptIterateAll(PyMOLGlobals *G)
{
  int active = false;
  float center_array[8] = {0.0F,0.0F,0.0F,0.0F};
  float *center = center_array;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  int state = SceneGetState(G);
  CGOReset(G->DebugCGO);

  if(SettingGet(G,cSetting_sculpting)) {
    if(!SettingGetGlobal_b(G,cSetting_sculpt_auto_center)) 
      center = NULL;
      
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptIterate(objMol,state,
                                      SettingGet_i(G,NULL,objMol->Obj.Setting,
                                                   cSetting_sculpting_cycles),
                                      center);
          active = true;
        }
      }
    }
    if(center && (center[3]>1.0F)) {
      float pos[3];
      SceneGetPos(G,pos);
      center[3] = 1.0F/center[3];
      scale3f(center,center[3],center);
      center[7] = 1.0F/center[7];
      scale3f(center+4,center[7],center+4);
      subtract3f(center,center+4,center);
      add3f(pos,center,center);
      ExecutiveCenter(G,NULL,-1,true,false,center,true);
    }
  }
  return(active);
}
/*========================================================================*/
float ExecutiveSculptIterate(PyMOLGlobals *G,char *name,int state,int n_cycle)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  register CExecutive *I = G->Executive;
  int ok=true;
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  float total_strain = 0.0F;

  if(state<0) state=SceneGetState(G);

  if(WordMatch(G,name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          total_strain+=ObjectMoleculeSculptIterate(objMol,state,n_cycle,NULL);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    total_strain=ObjectMoleculeSculptIterate((ObjectMolecule*)obj,state,n_cycle,NULL);
  }
  return(total_strain);
}
/*========================================================================*/
int ExecutiveSculptActivate(PyMOLGlobals *G,char *name,int state,int match_state,int match_by_segment)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  register CExecutive *I = G->Executive;
  int ok=true;
  if(state<0) state=SceneGetState(G);

  if(WordMatch(G,name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptImprint(objMol,state,match_state,match_by_segment);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    ObjectMoleculeSculptImprint((ObjectMolecule*)obj,state,match_state,match_by_segment);
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSculptDeactivate(PyMOLGlobals *G,char *name)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  SpecRec *rec = NULL;
  ObjectMolecule *objMol;
  register CExecutive *I = G->Executive;

  int ok=true;

  if(WordMatch(G,name,cKeywordAll,true)<0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          objMol =(ObjectMolecule*)rec->obj;
          ObjectMoleculeSculptClear(objMol);
        }
      }
    }
  } else if(!obj) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else if(obj->type != cObjectMolecule) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Executive-Error: object %s is not a molecular object.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    ObjectMoleculeSculptClear((ObjectMolecule*)obj);
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveSetGeometry(PyMOLGlobals *G,char *s1,int geom,int valence)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  int ok=false;

  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_SetGeometry;
    op1.i1 = geom;
    op1.i2 = valence;
    op1.i3 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    if(op1.i3) ok=true;
  } else {
    ErrMessage(G,"SetGeometry","Invalid selection.");
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveMultiSave(PyMOLGlobals *G,char *fname,char *name,int state,int append)
{
  int result=false;
  SpecRec *tRec;
  ObjectMolecule *objMol;
  
  PRINTFD(G,FB_Executive)
    " ExecutiveMultiSave-Debug: entered %s %s.\n",fname,name
    ENDFD;
  tRec = ExecutiveFindSpec(G,name);
  if(tRec) {
    if(tRec->type==cExecObject)
      if(tRec->obj->type==cObjectMolecule) {
        objMol =(ObjectMolecule*)tRec->obj;
        result = ObjectMoleculeMultiSave(objMol,fname,state,append);
      }
  }
  return(result);
  
}
int ExecutiveMapSetBorder(PyMOLGlobals *G,char *name,float level,int state)
{
  register CExecutive *I = G->Executive;
  int result=true;
  CTracker *I_Tracker= I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;

  while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
    if(rec) {
      switch(rec->type) {
      case cExecObject:
        if(rec->obj->type==cObjectMap) {
          ObjectMap *obj =(ObjectMap*)rec->obj;
          result = ObjectMapSetBorder(obj,level,state);
          
          if(result) { 
            ExecutiveInvalidateMapDependents(G,obj->Obj.Name);
          }
        }
        break;
      }
    }
  }
  
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveMapDouble(PyMOLGlobals *G,char *name,int state)
{
  register CExecutive *I = G->Executive;
  int result=true;
  CTracker *I_Tracker= I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;

  while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
    if(rec) {
      switch(rec->type) {
      case cExecObject:
        if(rec->obj->type==cObjectMap) {
          ObjectMap *obj =(ObjectMap*)rec->obj;
          result = ObjectMapDouble(obj,state);
          if(result) { 
            ExecutiveInvalidateMapDependents(G,obj->Obj.Name);
          }
          if(result && rec->visible)
            SceneChanged(G);
        }
        break;
      }
    }
  }

  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveMapHalve(PyMOLGlobals *G,char *name,int state,int smooth)
{
  register CExecutive *I = G->Executive;
  int result=true;
  CTracker *I_Tracker= I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  SpecRec *rec;

  while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
    if(rec) {
      switch(rec->type) {
      case cExecObject:
        if(rec->obj->type==cObjectMap) {
          ObjectMap *obj =(ObjectMap*)rec->obj;
          result = ObjectMapHalve(obj,state,smooth);
          if(result) { 
            ExecutiveInvalidateMapDependents(G,obj->Obj.Name);
          }
          if(result && rec->visible)
            SceneChanged(G);
        }
        break;
      }
    }
  }

  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
  return result;
}

int ExecutiveMapTrim(PyMOLGlobals *G,char *name,
                     char *sele,float buffer,
                     int map_state,int sele_state,
                     int quiet)
{
  register CExecutive *I = G->Executive;
  int result=true;
  float mn[3],mx[3];
  if(ExecutiveGetExtent(G,sele,mn,mx,true,sele_state,false)) {
    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    SpecRec *rec;

    {
      int a;
      float t;
      for(a=0;a<3;a++) {
        mn[a]-=buffer;
        mx[a]+=buffer;
        if(mn[a]>mx[a]) {
          t=mn[a]; mn[a]=mx[a]; mx[a]=t;
        }
      }
    }
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        switch(rec->type) {
        case cExecObject:
          if(rec->obj->type==cObjectMap) {
            ObjectMap *obj =(ObjectMap*)rec->obj;
            result = result && ObjectMapTrim(obj,map_state,mn,mx,quiet);
            if(result) 
              ExecutiveInvalidateMapDependents(G,obj->Obj.Name);
            if(result && rec->visible)
              SceneChanged(G);
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
  return result;
}

void ExecutiveSelectRect(PyMOLGlobals *G,BlockRect *rect,int mode)
{
  Multipick smp;
  OrthoLineType buffer,buf2;
  char selName[WordLength] = cLeftButSele;
  char prefix[3]="";
  int log_box = 0;
  int logging;
  char empty_string[1] = "";
  char *sel_mode_kw = empty_string;

  logging = (int)SettingGet(G,cSetting_logging);
  if(logging)
    log_box= (int)SettingGet(G,cSetting_log_box_selections);
  /*  if(logging==cPLog_pml)
      strcpy(prefix,"_ ");*/
  smp.picked=VLAlloc(Picking,1000);
  smp.x=rect->left;
  smp.y=rect->bottom;
  smp.w=rect->right-rect->left;
  smp.h=rect->top-rect->bottom;
  SceneMultipick(G,&smp);
  if(smp.picked[0].src.index) {
    SelectorCreate(G,cTempRectSele,NULL,NULL,1,&smp);
    if(log_box) SelectorLogSele(G,cTempRectSele);
    switch(mode) {
    case cButModeRect:
      if(mode==cButModeRect) {
        SelectorCreate(G,cLeftButSele,cTempRectSele,NULL,1,NULL);
        if(log_box) {
          sprintf(buf2,"%scmd.select(\"%s\",\"%s\",enable=1)\n",prefix,cLeftButSele,cTempRectSele);
          PLog(G,buf2,cPLog_no_flush);
        }
      } 
      break;
    case cButModeSeleAdd:
    case cButModeSeleSub:
      ExecutiveGetActiveSeleName(G,selName,true,SettingGet(G,cSetting_logging));
      sel_mode_kw = SceneGetSeleModeKeyword(G);        
      /* intentional omission of break! */
    case cButModeRectAdd:
    case cButModeRectSub:
      if(SelectorIndexByName(G,selName)>=0) {
        if((mode==cButModeRectAdd)||(mode==cButModeSeleAdd)) {
          sprintf(buffer,"(?%s or %s(%s))",selName,sel_mode_kw,cTempRectSele);
          SelectorCreate(G,selName,buffer,NULL,0,NULL);
          if(log_box) {
            sprintf(buf2,"%scmd.select(\"%s\",\"(%s)\",enable=1)\n",prefix,selName,buffer);
            PLog(G,buf2,cPLog_no_flush);
          }
        } else {
          sprintf(buffer,"(%s(?%s) and not %s(%s))",sel_mode_kw,selName,sel_mode_kw,cTempRectSele);
          SelectorCreate(G,selName,buffer,NULL,0,NULL);
          if(log_box) {
            sprintf(buf2,"%scmd.select(\"%s\",\"%s\",enable=1)\n",prefix,selName,buffer);
            PLog(G,buf2,cPLog_no_flush);
          }
        }
      } else {
        if((mode==cButModeRectAdd)||(mode==cButModeSeleAdd)) {
          sprintf(buffer,"%s(?%s)",sel_mode_kw,cTempRectSele);
          SelectorCreate(G,selName,buffer,NULL,0,NULL);
          if(log_box) {
            sprintf(buf2,"%scmd.select(\"%s\",\"%s\",enable=1)\n",prefix,selName,buffer);
            PLog(G,buf2,cPLog_no_flush);
          }
        } else {
          SelectorCreate(G,selName,"(none)",NULL,0,NULL);
          if(log_box) {
            sprintf(buf2,"%scmd.select(\"%s\",\"(none)\",enable=1)\n",prefix,selName);
            PLog(G,buf2,cPLog_no_flush);
          }
        }
      }
      if(SettingGet(G,cSetting_auto_show_selections)) {
        ExecutiveSetObjVisib(G,selName,true,false);
      }
      break;
    }
    if(log_box) {
      sprintf(buf2,"%scmd.delete(\"%s\")\n",prefix,cTempRectSele);
      PLog(G,buf2,cPLog_no_flush);
      PLogFlush(G);
    }
    ExecutiveDelete(G,cTempRectSele);
    WizardDoSelect(G,selName);
  }
  VLAFreeP(smp.picked);
}

int ExecutiveTranslateAtom(PyMOLGlobals *G,char *sele,float *v,int state,int mode,int log)
{
  int ok=true;
  ObjectMolecule *obj0;
  int sele0 = SelectorIndexByName(G,sele);
  int i0;
  if(sele0<0) {
    PRINTFB(G,FB_Executive,FB_Errors)
      "Error: bad selection %s.\n",sele
      ENDFB(G);
    ok=false;
  } else{ 
    obj0 = SelectorGetSingleObjectMolecule(G,sele0);
    if(!obj0) {
      PRINTFB(G,FB_Executive,FB_Errors)
        "Error: selection isn't a single atom.\n"
        ENDFB(G);
      ok=false;
    } else {
      i0 = ObjectMoleculeGetAtomIndex(obj0,sele0);
      if(i0<0) {
        PRINTFB(G,FB_Executive,FB_Errors)
          "Error: selection isn't a single atom.\n"
          ENDFB(G);
        ok=false;
      } else {
        ObjectMoleculeMoveAtom(obj0,state,i0,v,mode,log);
      }
    }
  }
  return(ok);
}

int ExecutiveCombineObjectTTT(PyMOLGlobals *G,char *name,float *ttt, int reverse_order)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  int ok=true;

  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    ObjectCombineTTT(obj,ttt,reverse_order);
    if(obj->fInvalidate)
      obj->fInvalidate(obj,cRepNone,cRepInvExtents,-1);
    SceneInvalidate(G);
  }
  return(ok);
}

int ExecutiveSetObjectTTT(PyMOLGlobals *G,char *name,float *ttt,int state,int quiet)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  int ok=true;

  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    ObjectSetTTT(obj,ttt,state);
    if(obj->fInvalidate)
      obj->fInvalidate(obj,cRepNone,cRepInvExtents,-1);
  }
  return(ok);
}

int ExecutiveGetObjectTTT(PyMOLGlobals *G,char *name,float **ttt,int state,int quiet)
{
  CObject *obj = ExecutiveFindObjectByName(G,name);
  int ok=true;

  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
    ok=false;
  } else {
    ObjectGetTTT(obj,ttt,state);
  }
  return(ok);
}

int ExecutiveTransformSelection(PyMOLGlobals *G,int state,char *s1,int log,float *ttt,int homogenous)
{
  int sele=-1;
  ObjectMolecule *obj = NULL;
  ObjectMolecule **vla = NULL;
  int nObj;
  int ok=true;
  int a;

  sele = SelectorIndexByName(G,s1);
  if(sele<0)
    ok=false;
  if(ok) {
    vla=SelectorGetObjectMoleculeVLA(G,sele);
    if(!vla) ok=false;
  }
  if(ok) {
    nObj = VLAGetSize(vla);
    for(a=0;a<nObj;a++) {
      obj=vla[a];
      ObjectMoleculeTransformSelection(obj,state,sele,ttt,log,s1,homogenous,true);
    }
  }
  SceneInvalidate(G);
  VLAFreeP(vla);
  return(ok);
}

int ExecutiveTransformObjectSelection2(PyMOLGlobals *G,CObject *obj,int state,
                                       char *s1,int log,float *matrix,
                                       int homogenous,int global)
{
  int ok=true;
  
  switch(obj->type) {
  case cObjectMolecule: 
    {
      int sele=-1;
      ObjectMolecule *objMol = (ObjectMolecule*)obj;
      
      if(s1&&s1[0]) {
        sele = SelectorIndexByName(G,s1);
        if(sele<0)
          ok=false;
      }
      if(!ok) {
        PRINTFB(G,FB_ObjectMolecule,FB_Errors)
          "Error: selection object %s not found.\n",s1
          ENDFB(G);
      } else {
        ObjectMoleculeTransformSelection(objMol,state,sele,matrix,log,s1,homogenous,global);
      }
      EditorDihedralInvalid(G,objMol);
      SceneInvalidate(G);
    }
    break;
  case cObjectMap:
    {
      double matrixd[116];
      if(homogenous) {
        convert44f44d(matrix,matrixd);
      } else {
        convertTTTfR44d(matrix,matrixd);
      }
      ObjectMapTransformMatrix((ObjectMap*)obj,state,matrixd);
    }
    break;   
  case cObjectGroup:
    {
      double matrixd[116];
      if(homogenous) {
        convert44f44d(matrix,matrixd);
      } else {
        convertTTTfR44d(matrix,matrixd);
      }
      ObjectGroupTransformMatrix((ObjectGroup*)obj,state,matrixd);
    }
    break;   
  }
  return(ok);
}

int ExecutiveTransformObjectSelection(PyMOLGlobals *G,char *name,int state,
                                      char *s1,int log,float *matrix,
                                      int homogenous,int global)
{
  int ok=true;

  CObject *obj = ExecutiveFindObjectByName(G,name);
  if(obj) {
    return ExecutiveTransformObjectSelection2(G,obj,state,s1,log,matrix,homogenous,global);
  }
  return ok;
}

int ExecutiveValidName(PyMOLGlobals *G,char *name)
{
  int result=true;
  if(!ExecutiveFindSpec(G,name)) {
    int ignore_case = SettingGetGlobal_b(G,cSetting_ignore_case);

    if(!WordMatchExact(G,name,cKeywordAll,ignore_case))
      if(!WordMatchExact(G,name,cKeywordSame,ignore_case))
        if(!WordMatchExact(G,name,cKeywordCenter,ignore_case))
          if(!WordMatchExact(G,name,cKeywordOrigin,ignore_case))
            result=false;
  }
  return result;
}

int ExecutivePhiPsi(PyMOLGlobals *G,char *s1,ObjectMolecule ***objVLA,int **iVLA,
                    float **phiVLA,float **psiVLA,int state) 
{
  int sele1=SelectorIndexByName(G,s1);
  int result = false;
  ObjectMoleculeOpRec op1;
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.i1 = 0;
    op1.i2 = state;
    op1.obj1VLA=VLAlloc(ObjectMolecule*,1000);
    op1.i1VLA=VLAlloc(int,1000);
    op1.f1VLA=VLAlloc(float,1000);
    op1.f2VLA=VLAlloc(float,1000);
    op1.code=OMOP_PhiPsi;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    result = op1.i1;
    VLASize(op1.i1VLA,int,op1.i1);
    VLASize(op1.obj1VLA,ObjectMolecule*,op1.i1);
    VLASize(op1.f1VLA,float,op1.i1);
    VLASize(op1.f2VLA,float,op1.i1);
    *iVLA=op1.i1VLA;
    *objVLA=op1.obj1VLA;
    *phiVLA=op1.f1VLA;
    *psiVLA=op1.f2VLA;
  } else {
    *objVLA=NULL;
    *iVLA=NULL;
    *phiVLA=NULL;
    *psiVLA=NULL;
  }
  return(result);
}


int ExecutiveAlign(PyMOLGlobals *G,char *s1,char *s2,char *mat_file,float gap, float extend,
                   int max_gap, int max_skip, float cutoff,int cycles,int quiet,char *oname,
                   int state1,int state2, ExecutiveRMSInfo *rms_info,int transform,int reset,
                   float seq_wt,float radius, float scale, float base, float coord_wt,
                   float expect, int window, float ante)
{
  int sele1=SelectorIndexByName(G,s1);
  int sele2=SelectorIndexByName(G,s2);
  int *vla1=NULL;
  int *vla2=NULL;
  int na,nb;
  int c;
  int ok=true;
  int use_sequence = (mat_file && mat_file[0] && (seq_wt!=0.0F));
  int use_structure = (seq_wt>=0.0F); /* negative seq_wt means sequence only! */
  ObjectMolecule *mobile_obj = NULL;
  CMatch *match = NULL;

  if(!use_structure) window = 0;
  
  if((scale==0.0F) && (seq_wt==0.0F) && (ante<0.0F) && window)
    ante = window;

  if(ante<0.0F)
    ante=0.0F;

  if((sele1>=0)) {
    mobile_obj = SelectorGetSingleObjectMolecule(G,sele1);
    if(!mobile_obj) {
      ok=false;
      PRINTFB(G,FB_Executive,FB_Errors)
        " ExecutiveAlign: mobile selection must derive from one object only.\n"
        ENDFB(G);
    }
  }
  if(ok&&(sele1>=0)&&(sele2>=0)&&rms_info) {
    vla1=SelectorGetResidueVLA(G,sele1,use_structure,NULL);
    vla2=SelectorGetResidueVLA(G,sele2,use_structure,mobile_obj);
    if(vla1&&vla2) {
      na = VLAGetSize(vla1)/3;
      nb = VLAGetSize(vla2)/3;
      if(na&&nb) {
        match = MatchNew(G,na,nb,window);
        if(match) {
          if(use_sequence) {
            if (ok) ok = MatchResidueToCode(match,vla1,na);
            if (ok) ok = MatchResidueToCode(match,vla2,nb);
            if (ok) ok = MatchMatrixFromFile(match,mat_file,quiet);
            if (ok) ok = MatchPreScore(match,vla1,na,vla2,nb,quiet);
          }
          if(use_structure) {
            if (ok) ok = SelectorResidueVLAsTo3DMatchScores(G,match,
                                                            vla1,na,state1,
                                                            vla2,nb,state2,seq_wt,
                                                            radius,scale,base,
                                                            coord_wt,expect);
          }
          if(ok) ok = MatchAlign(match,gap,extend,max_gap,
                                 max_skip,quiet,window,ante);
          if(ok) {
            rms_info->raw_alignment_score = match->score;
            rms_info->n_residues_aligned = match->n_pair;
            if(match->pair) { 
            
              c = SelectorCreateAlignments(G,match->pair,
                                           sele1,vla1,sele2,vla2,
                                           "_align1","_align2",false,false);

              if(c) {
                int mode = 2;
                if(!quiet) {
                  PRINTFB(G,FB_Executive,FB_Actions)
                    " ExecutiveAlign: %d atoms aligned.\n",c
                    ENDFB(G);
                }
                if(oname&&oname[0]&&reset)
                  ExecutiveDelete(G,oname);
                if(!transform)
                  mode = 1;
                ok = ExecutiveRMS(G,"_align1","_align2",mode,cutoff,cycles,
                                  quiet,oname,
                                  state1,state2,false,0, rms_info);
              }
            }
          }
          MatchFree(match);
        }
      } else {
        ok=false;
        PRINTFB(G,FB_Executive,FB_Errors)
          " ExecutiveAlign: invalid selections for alignment.\n"
          ENDFB(G);
      }
    }
  }
  
  VLAFreeP(vla1);
  VLAFreeP(vla2);
  return ok;
}

int ExecutivePairIndices(PyMOLGlobals *G,char *s1,char *s2,int state1,int state2,
                         int mode,float cutoff,float h_angle,
                         int **indexVLA, ObjectMolecule ***objVLA)
{
  int result = 0;
  int sele1,sele2;

  sele1 = SelectorIndexByName(G,s1);
  sele2 = SelectorIndexByName(G,s2);
  if((sele1>=0)&&(sele2>=0)) {
    result=SelectorGetPairIndices(G,sele1,state1,sele2,state2,
                                  mode,cutoff,h_angle,indexVLA,objVLA);
  } else {
    ErrMessage(G,"ExecutivePairIndices","One or more bad selections.");
  }
  return(result);
}


int ExecutiveCartoon(PyMOLGlobals *G,int type,char *s1)
{
  int sele1;
  ObjectMoleculeOpRec op1;

  sele1=SelectorIndexByName(G,s1);
  ObjectMoleculeOpRecInit(&op1);
  op1.i2=0;
  if(sele1>=0) {
    op1.code=OMOP_INVA;
    op1.i1=cRepCartoon; 
    op1.i2=cRepInvRep;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    op1.code = OMOP_Cartoon;
    op1.i1 = type;
    op1.i2 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
  } else {
    ErrMessage(G,"Cartoon","Invalid selection.");
  }
  return(op1.i2);
}
/*========================================================================*/
float *ExecutiveGetVertexVLA(PyMOLGlobals *G,char *s1,int state)
{
  /* returns NULL if none found */

  float *result = NULL;
  ObjectMoleculeOpRec op1;
  int sele1;
  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.nvv1 = 0;
    op1.vv1=VLAlloc(float,1000);
    if(state>=0) {
      op1.cs1 = state;
      op1.code=OMOP_SingleStateVertices;
    } else {
      op1.code=OMOP_VERT;
    }
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    VLASize(op1.vv1,float,op1.nvv1*3);
    result = op1.vv1;
  }
  return(result);
}
/*========================================================================*/
PyObject *ExecutiveGetSettingOfType(PyMOLGlobals *G,int index,
                                    char *object,int state,int type)
{ 
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* Assumes blocked Python interpreter */
  PyObject *result = NULL;
  CObject *obj = NULL;
  CSetting **handle=NULL,*set_ptr1=NULL,*set_ptr2=NULL;
  int ok=true;

  if(object)
    if(object[0]) {
      obj=ExecutiveFindObjectByName(G,object);
      if(!obj) 
        ok=false;
    } 
  if(!ok) {
    PRINTFB(G,FB_Executive,FB_Errors)
      " SettingGet-Error: object \"%s\" not found.\n",object
      ENDFB(G);
    ok=false;
  } else if(obj) {
    handle = obj->fGetSettingHandle(obj,-1);
    if(handle) set_ptr1 = *handle;
    if(state>=0) {
      handle = obj->fGetSettingHandle(obj,state);
      if(handle) 
        set_ptr2 = *handle;
      else {
        PRINTFB(G,FB_Executive,FB_Errors)
          " SettingGet-Error: object \"%s\" lacks state %d.\n",object,state+1
          ENDFB(G);
        ok=false;
      }
    }
  }
  if(ok) {
    switch(type) {
    case cSetting_boolean:
      {
        int value = SettingGet_b(G,set_ptr2,set_ptr1,index);
        result=Py_BuildValue("i",value);
      }
      break;
    case cSetting_int:
      {
        int value = SettingGet_i(G,set_ptr2,set_ptr1,index);
        result=Py_BuildValue("i",value);
      }
      break;
    case cSetting_float:
      {
        float value = SettingGet_f(G,set_ptr2,set_ptr1,index);
        result=Py_BuildValue("f",value);
      }
      break;
    case cSetting_float3:
      {
        float value[3];
        SettingGet_3f(G,set_ptr2,set_ptr1,index,value);
        result=Py_BuildValue("fff",value[0],value[1],value[2]);
      }
      break;
    case cSetting_color:
      {
        int value = SettingGet_color(G,set_ptr2,set_ptr1,index);
        result=Py_BuildValue("i",value);
      }
      break;
    case cSetting_string:
      {  
        OrthoLineType buffer = "";
        buffer[0]=0;
        SettingGetTextValue(G,set_ptr2,set_ptr1,index,buffer);
        result=Py_BuildValue("s",buffer);
      }
      break;
    default:
      result=Py_BuildValue("i",0);
      break;
    }
  } 
  return(result);
#endif

}
/*========================================================================*/
PyObject *ExecutiveGetSettingText(PyMOLGlobals *G,int index,char *object,int state)
{ 
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* Assumes blocked Python interpreter */
  PyObject *result = NULL;
  OrthoLineType buffer = "";
  CObject *obj = NULL;
  CSetting **handle=NULL,*set_ptr1=NULL,*set_ptr2=NULL;
  int ok=true;

  if(object)
    if(object[0]) {
      obj=ExecutiveFindObjectByName(G,object);
      if(!obj) 
        ok=false;
    } 
  if(!ok) {
    PRINTFB(G,FB_Executive,FB_Errors)
      " SettingGet-Error: object \"%s\" not found.\n",object
      ENDFB(G);
    ok=false;
  } else if(obj) {
    handle = obj->fGetSettingHandle(obj,-1);
    if(handle) set_ptr1 = *handle;
    if(state>=0) {
      handle = obj->fGetSettingHandle(obj,state);
      if(handle) 
        set_ptr2 = *handle;
      else {
        PRINTFB(G,FB_Executive,FB_Errors)
          " SettingGet-Error: object \"%s\" lacks state %d.\n",object,state+1
          ENDFB(G);
        ok=false;
      }
    }
  }
  if(ok) {
    buffer[0]=0;
    SettingGetTextValue(G,set_ptr2,set_ptr1,index,buffer);
    result=Py_BuildValue("s",buffer);
  } 
  
  return(result);
#endif

}
/*========================================================================*/
PyObject *ExecutiveGetSettingTuple(PyMOLGlobals *G,int index,char *object,int state)
{ /* Assumes blocked Python interpreter */
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;
  CSetting **handle = NULL;
  CObject *obj=NULL;
  int ok = true;
  PRINTFD(G,FB_Executive)
    " ExecutiveGetSettingTuple: object %p state %d\n",object,state
    ENDFD;

  if(object[0]==0) /* global */
    result = SettingGetTuple(G,NULL,NULL,index);
  else {

    if(strlen(object)) {
      obj=ExecutiveFindObjectByName(G,object);
      if(!obj) 
        ok=false;
    } else ok=false;
    if(!ok) {
      PRINTFB(G,FB_Executive,FB_Errors)
        " Executive: object not found.\n"
        ENDFB(G);
    } else {
      handle = obj->fGetSettingHandle(obj,state);
      if(handle) 
        result = SettingGetDefinedTuple(G,*handle,index);      
    }
  }
  if(!ok) {
    result = PConvAutoNone(Py_None);
  }
  return(result);
#endif
}
/*========================================================================*/
void ExecutiveSetLastObjectEdited(PyMOLGlobals *G,CObject *o)
{
  register CExecutive *I = G->Executive;
  I->LastEdited = o;
}
/*========================================================================*/
CObject *ExecutiveGetLastObjectEdited(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  return(I->LastEdited);
}
/*========================================================================*/
int ExecutiveSaveUndo(PyMOLGlobals *G,char *s1,int state)
{
  int sele1;
  ObjectMoleculeOpRec op1;

  if(state<0) state = SceneGetState(G);                
  sele1=SelectorIndexByName(G,s1);
  ObjectMoleculeOpRecInit(&op1);
  op1.i2=0;
  if(sele1>=0) {
    op1.code = OMOP_SaveUndo;
    op1.i1 = state;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
  }
  return(op1.i2);
}

/*========================================================================*/
int ExecutiveSetTitle(PyMOLGlobals *G,char *name,int state,char *text)
{
  int result=false;
  ObjectMolecule *obj;
  obj =ExecutiveFindObjectMoleculeByName(G,name);
  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
  } else {
    result = ObjectMoleculeSetStateTitle(obj,state,text);
  }
  SceneDirty(G);
  return(result);
}
/*========================================================================*/
char *ExecutiveGetTitle(PyMOLGlobals *G,char *name,int state)
{
  char *result = NULL;
  ObjectMolecule *obj;
  obj =ExecutiveFindObjectMoleculeByName(G,name);
  if(!obj) {
    PRINTFB(G,FB_ObjectMolecule,FB_Errors)
      "Error: object %s not found.\n",name 
      ENDFB(G);
  } else {
    result = ObjectMoleculeGetStateTitle(obj,state);
  }
  return(result);
}
/*========================================================================*/
void ExecutiveHideSelections(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection) {
      if(rec->visible) {
        rec->visible=false;
        SceneInvalidate(G);
        SeqDirty(G);
      }
    }
  }
}
/*========================================================================*/
void ExecutiveRenderSelections(PyMOLGlobals *G,int curState)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  int any_active = false;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecSelection) {
      if(rec->visible) {
        any_active = true;
        break;
      }
    }
  }

  if(any_active) {
    SpecRec *rec1;
    int sele;
    int no_depth;
    float min_width;
    float gl_width;
    int width;
    
    int max_width = (int)SettingGetGlobal_f(G,cSetting_selection_width_max);
    float width_scale = SettingGetGlobal_f(G,cSetting_selection_width_scale);
    int round_points = SettingGetGlobal_b(G,cSetting_selection_round_points);
    int vis_only = SettingGetGlobal_b(G,cSetting_selection_visible_only);
    int fog = SettingGet(G,cSetting_depth_cue)&&SettingGet(G,cSetting_fog);
    rec = NULL;
    min_width = SettingGetGlobal_f(G,cSetting_selection_width);
    
    if(width_scale>=0.0F) {
      width = (int)((width_scale*SettingGetGlobal_f(G,cSetting_stick_radius)/
                     SceneGetScreenVertexScale(G,NULL)));
      if(width<min_width)
        width = (int)min_width;
      else if(width>max_width)
        width = (int)max_width;
    } else
      width = (int)min_width;
    
    if(round_points) {
      glEnable(GL_POINT_SMOOTH);
      glAlphaFunc(GL_GREATER, 0.5F);
      glEnable(GL_ALPHA_TEST);
      glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
      width = (int)(width*1.44F);
    } else {
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_ALPHA_TEST);
      glHint(GL_POINT_SMOOTH_HINT,GL_FASTEST);
    }
    
    no_depth = (int)SettingGet(G,cSetting_selection_overlay);
    
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecSelection) {
        
        if(rec->visible) {
          int enabled = true;
          SpecRec *group_rec = rec->group;
          while(enabled && group_rec ) { 
            if(!group_rec->visible)
              enabled=false;
            else
              group_rec = group_rec->group;
          }
          
          if(enabled) {
            
            sele = SelectorIndexByName(G,rec->name); /* TODO: speed this up */
            if(sele>=0) {
              
              if(no_depth)
                glDisable(GL_DEPTH_TEST);
              glDisable(GL_FOG);
              
              if(rec->sele_color<0)
                glColor3f(1.0F,0.2F,0.6F);
              else
                glColor3fv(ColorGet(G,rec->sele_color));
            
              gl_width=(float)width;
              if(width>6) { /* keep it even above 6 */
                if(width&0x1) {
                  width--;
                  gl_width = (float)width;
                }
              }
              glPointSize(gl_width);
              glBegin(GL_POINTS);
              rec1 = NULL;
              while(ListIterate(I->Spec,rec1,next)) {
                if(rec1->type==cExecObject) {
                  if(rec1->obj->type==cObjectMolecule) {
                    ObjectMoleculeRenderSele((ObjectMolecule*)rec1->obj,curState,sele,vis_only);
                  }
                }
              }
              glEnd();
            
              if(width>2) {
                switch(width) {
                case 1:
                case 2:
                case 3:
                  glPointSize(1.0F);
                  break;
                case 4:
                  glPointSize(2.0F);
                  break;
                case 5:
                  glPointSize(3.0F); 
                  break;
                case 6:
                case 7:
                case 8:
                case 9:
                  glPointSize(4.0F);
                  break;
                default:
                  glPointSize(6.0F);
                  break;
                }
              
                glColor3f(0.0F,0.0F,0.0F);
                glBegin(GL_POINTS);
                rec1 = NULL;
                while(ListIterate(I->Spec,rec1,next)) {
                  if(rec1->type==cExecObject) {
                    if(rec1->obj->type==cObjectMolecule) {
                      ObjectMoleculeRenderSele((ObjectMolecule*)rec1->obj,curState,sele,vis_only);
                    }
                  }
                }
                glEnd();
              }
            
              if(width>4) {
                if(width>5) {
                  glPointSize(2.0F);
                }
                else 
                  glPointSize(1.0F);
                glColor3f(1.0F,1.0F,1.0F);
              
                glBegin(GL_POINTS);
                rec1 = NULL;
                while(ListIterate(I->Spec,rec1,next)) {
                  if(rec1->type==cExecObject) {
                    if(rec1->obj->type==cObjectMolecule) {
                      ObjectMoleculeRenderSele((ObjectMolecule*)rec1->obj,curState,sele,vis_only);
                    }
                  }
                }
                glEnd();
              }
            
            
            
              if(no_depth)
                glEnable(GL_DEPTH_TEST);
              if(fog)
                glEnable(GL_FOG);
            }
          }
        }
      }
    }
    if(round_points) {
      glAlphaFunc(GL_GREATER, 0.05F);
    }
  }
}
/*========================================================================*/
int ExecutiveGetDistance(PyMOLGlobals *G,char *s0,char *s1,float *value,int state)
{
  /* TO DO: add support for averaging over multiple states */

  Vector3f v0,v1;
  int sele0=-1,sele1=-1;
  int ok=true;

  if((sele0 = SelectorIndexByName(G,s0))<0)
    ok = ErrMessage(G,"GetDistance","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(G,s1))<0)
    ok = ErrMessage(G,"GetDistance","Selection 2 invalid.");    
  if(ok) {
    if (!SelectorGetSingleAtomVertex(G,sele0,state,v0))
      ok = ErrMessage(G,"GetDistance","Selection 1 doesn't contain a single atom/vertex.");
    if (!SelectorGetSingleAtomVertex(G,sele1,state,v1))
      ok = ErrMessage(G,"GetDistance","Selection 2 doesn't contain a single atom/vertex.");
  }
  if(ok) {
    (*value)=(float)diff3f(v0,v1);
  }
  return ok;
}
/*========================================================================*/
int ExecutiveGetAngle(PyMOLGlobals *G,char *s0,char *s1,char *s2,float *value,int state)
{

  /* TO DO: add support for averaging over multiple states */

  Vector3f v0,v1,v2;
  int sele0=-1,sele1=-1,sele2=-1;
  int ok=true;
  float d1[3],d2[3];
  if((sele0 = SelectorIndexByName(G,s0))<0)
    ok = ErrMessage(G,"GetAngle","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(G,s1))<0)
    ok = ErrMessage(G,"GetAngle","Selection 2 invalid.");    
  else if((sele2 = SelectorIndexByName(G,s2))<0)
    ok = ErrMessage(G,"GetAngle","Selection 3 invalid.");
  if(ok) {
    if (!SelectorGetSingleAtomVertex(G,sele0,state,v0))
      ok = ErrMessage(G,"GetAngle","Selection 1 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele1,state,v1))
      ok = ErrMessage(G,"GetAngle","Selection 2 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele2,state,v2))
      ok = ErrMessage(G,"GetAngle","Selection 3 doesn't contain a single atom/vertex.");          
  }
  if(ok) {
    subtract3f(v0,v1,d1);
    subtract3f(v2,v1,d2);
    (*value)=rad_to_deg(get_angle3f(d1,d2));
  }
  return ok;
}
/*========================================================================*/
int ExecutiveGetDihe(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,float *value,int state)
{

  /* TO DO: add support for averaging over multiple states */

  Vector3f v0,v1,v2,v3;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int ok=true;
  
  if((sele0 = SelectorIndexByName(G,s0))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(G,s1))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 2 invalid.");    
  else if((sele2 = SelectorIndexByName(G,s2))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 3 invalid.");
  else if((sele3 = SelectorIndexByName(G,s3))<0)
    ok = ErrMessage(G,"GetDihedral","Selection 4 invalid.");
  if(ok) {
    if (!SelectorGetSingleAtomVertex(G,sele0,state,v0))
      ok = ErrMessage(G,"GetDihedral","Selection 1 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele1,state,v1))
      ok = ErrMessage(G,"GetDihedral","Selection 2 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele2,state,v2))
      ok = ErrMessage(G,"GetDihedral","Selection 3 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele3,state,v3))
      ok = ErrMessage(G,"GetDihedral","Selection 4 doesn't contain a single atom/vertex.");          
  }
  if(ok) {
    (*value)=rad_to_deg(get_dihedral3f(v0,v1,v2,v3));
  }
  return ok;
}
/*========================================================================*/
int ExecutiveSetDihe(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,float value,int state,int quiet)
{
  Vector3f v0,v1,v2,v3;
  int sele0=-1,sele1=-1,sele2=-1,sele3=-1;
  int ok=true;
  int save_state;
  float current;
  float change;

  if((sele0 = SelectorIndexByName(G,s0))<0)
    ok = ErrMessage(G,"SetDihedral","Selection 1 invalid.");    
  else if((sele1 = SelectorIndexByName(G,s1))<0)
    ok = ErrMessage(G,"SetDihedral","Selection 2 invalid.");    
  else if((sele2 = SelectorIndexByName(G,s2))<0)
    ok = ErrMessage(G,"SetDihedral","Selection 3 invalid.");
  else if((sele3 = SelectorIndexByName(G,s3))<0)
    ok = ErrMessage(G,"SetDihedral","Selection 4 invalid.");
  if(ok) {
    if (!SelectorGetSingleAtomVertex(G,sele0,state,v0))
      ok = ErrMessage(G,"SetDihedral","Selection 1 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele1,state,v1))
      ok = ErrMessage(G,"SetDihedral","Selection 2 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele2,state,v2))
      ok = ErrMessage(G,"SetDihedral","Selection 3 doesn't contain a single atom/vertex.");          
    if (!SelectorGetSingleAtomVertex(G,sele3,state,v3))
      ok = ErrMessage(G,"SetDihedral","Selection 4 doesn't contain a single atom/vertex.");          
  }
  if(ok) {
    current=rad_to_deg(get_dihedral3f(v0,v1,v2,v3));
    change=value-current;
    save_state = SceneGetState(G);                
    SceneSetFrame(G,-1,state); /* KLUDGE ALERT!
                                * necessary because the editor 
                                * can only work on the current state...this
                                * needs to be changed.*/
    EditorSelect(G,s2,s1,NULL,NULL,false,true,true);
    EditorTorsion(G,change);
    SceneSetFrame(G,-1,save_state);
    if(!quiet) {
      PRINTFB(G,FB_Editor,FB_Actions)
        " SetDihedral: adjusted to %5.3f\n",value
        ENDFB(G);
    }

  }
  return ok;
}
/*========================================================================*/
float ExecutiveGetArea(PyMOLGlobals *G,char *s0,int sta0,int load_b)
{
  ObjectMolecule *obj0;
  RepDot *rep;
  CoordSet *cs;
  float result=-1.0F;
  int a,sele0;
  int known_member=-1;
  int is_member;
  int *ati;
  float *area;
  AtomInfoType *ai=NULL;
  ObjectMoleculeOpRec op;
  sele0 = SelectorIndexByName(G,s0);
  if(sele0<0) {
    ErrMessage(G,"Area","Invalid selection.");
  } else {
    obj0 = SelectorGetSingleObjectMolecule(G,sele0);
    if(!(obj0)) {
      if(SelectorCountAtoms(G,sele0,sta0)>0)
        ErrMessage(G,"Area","Selection must be within a single object.");
      else
        result = 0.0F;
    }
    else {
      cs = ObjectMoleculeGetCoordSet(obj0,sta0);
      if(!cs)
        ErrMessage(G,"Area","Invalid state.");
      else {
        rep = (RepDot*)RepDotDoNew(cs,cRepDotAreaType,sta0);
        if(!rep) 
          ErrMessage(G,"Area","Can't get dot representation.");
        else {

          if(load_b) {
            /* zero out B-values within selection */
            ObjectMoleculeOpRecInit(&op);
            op.code=OMOP_SetB;
            op.f1=0.0;
            op.i1=0;
            ExecutiveObjMolSeleOp(G,sele0,&op);
          }

          result=0.0;
          
          area=rep->A;
          ati=rep->Atom;
          
          is_member = false;

          for(a=0;a<rep->N;a++) {
            
            if(known_member!=(*ati)) {
              known_member=(*ati);
              ai=obj0->AtomInfo+known_member;
              is_member = SelectorIsMember(G,ai->selEntry,sele0);
            } 

            if(is_member) {
              result+=(*area);
              if(load_b)
                ai->b+=(*area);
            }
            area++;
            ati++;
          }
          
          rep->R.fFree((Rep*)rep); /* free the representation */
        }
      }
    }
  }
  return(result);
}

/*========================================================================*/
char *ExecutiveGetNames(PyMOLGlobals *G,int mode,int enabled_only,char *s0)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  char *result;
  int size = 0;
  int stlen;
  int sele0 = -1;
  int incl_flag = 0;
  if(s0[0]) {
    sele0 = SelectorIndexByName(G,s0);    
  }
  result=VLAlloc(char,1000);

  while(ListIterate(I->Spec,rec,next)) {
    if(
       (rec->type==cExecObject&&((!mode)||(mode==1)||(mode==3)||(mode==4)))||
       (rec->type==cExecSelection&&((!mode)||(mode==2)||(mode==3)||(mode==5))))
      {
        if((mode<3)||(rec->name[0]!='_')) {
          if((!enabled_only)||(rec->visible)) {
            stlen = strlen(rec->name);
            if(sele0<0) 
              incl_flag = 1;
            else 
              switch(rec->type) {
              case cExecObject:
                if(rec->obj->type == cObjectMolecule) {
                  int a;
                  ObjectMolecule *obj_mol = (ObjectMolecule*)rec->obj;
                  AtomInfoType *ai = obj_mol->AtomInfo;
                  for(a=0;a<obj_mol->NAtom;a++) {
                    if(SelectorIsMember(G,ai->selEntry,sele0)) {
                      incl_flag = 1;
                      break;
                    }
                    ai++;
                  }
                }
                break;
              case cExecSelection:
                if(SelectorCheckIntersection(G,sele0,SelectorIndexByName(G,rec->name))) {
                  incl_flag=1;
                  break;
                }
                break;
              }
            if(incl_flag) {
              VLACheck(result,char,size+stlen+1);
              strcpy(result+size,rec->name);
              size+=stlen+1;
            }
          }
        }
      }
    
  }
  VLASize(result,char,size);
  return(result);
}
/*========================================================================*/
int ExecutiveGetType(PyMOLGlobals *G,char *name,WordType type)
{
  SpecRec *rec = NULL;
  int ok=true;
  rec = ExecutiveFindSpec(G,name);
  if(!rec) {
    ok=false;
  } else {
    if(rec->type==cExecObject) {
      strcpy(type,"object:");
      if(rec->obj->type==cObjectMolecule)
        strcat(type,"molecule");
      else if(rec->obj->type==cObjectMap)
        strcat(type,"map");
      else if(rec->obj->type==cObjectMesh)
        strcat(type,"mesh");
      else if(rec->obj->type==cObjectSlice)
        strcat(type,"slice");
      else if(rec->obj->type==cObjectSurface)
        strcat(type,"surface");
      else if(rec->obj->type==cObjectMeasurement)
        strcat(type,"measurement");
      else if(rec->obj->type==cObjectCGO)
        strcat(type,"cgo");
      else if(rec->obj->type==cObjectGroup)
        strcat(type,"group");
    } else if(rec->type==cExecSelection) {
      strcpy(type,"selection");
    }
  }
  return(ok);
}

/*========================================================================*/
void ExecutiveUpdateCmd(PyMOLGlobals *G,char *s0,char *s1,int sta0,int sta1,
                        int method, int quiet)
{
  int sele0,sele1;

  sele0 = SelectorIndexByName(G,s0);
  sele1 = SelectorIndexByName(G,s1);
  if(!(sele0&&sele1)) {
    ErrMessage(G,"Update","One or more invalid input selections.");
  } else {
    SelectorUpdateCmd(G,sele0,sele1,sta0,sta1,method,quiet);
  }
}
/*========================================================================*/
void ExecutiveRenameObjectAtoms(PyMOLGlobals *G,char *name,int force) 
{
  register CExecutive *I = G->Executive;
  CObject *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(G,name);
    if(!os)
      ErrMessage(G," Executive","object not found.");
    else if(os->type!=cObjectMolecule) {
      ErrMessage(G," Executive","bad object type.");
      os = NULL;
    }
  }
  
  if(os||(!strlen(name))) { /* sort one or all */
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule)
          if((!os)||(rec->obj==os)) {
            obj =(ObjectMolecule*)rec->obj;
            ObjectMoleculeRenameAtoms(obj,force);  
          }
    }
    SceneChanged(G);
  }
} 

/*========================================================================*/
int  ExecutiveInvert(PyMOLGlobals *G,int quiet)
{
  int ok=false;
  ok = EditorInvert(G,quiet);
  return(ok);
}
/*========================================================================*/
void ExecutiveFuse(PyMOLGlobals *G,char *s0,char *s1,int mode,int recolor,int move_flag)
{
  int i0=-1;
  int i1=-1;
  int sele0,sele1,sele2;
  ObjectMolecule *obj0,*obj1;
  ObjectMoleculeOpRec op;
  
#define tmp_fuse_sele "tmp_fuse_sele"

  sele0 = SelectorIndexByName(G,s0);
  if(sele0>=0) {
    sele1 = SelectorIndexByName(G,s1);
    if(sele1>=0) {
      EditorInactivate(G);
      obj0 = SelectorGetSingleObjectMolecule(G,sele0);
      obj1 = SelectorGetSingleObjectMolecule(G,sele1);
      if(obj0)
        i0 = ObjectMoleculeGetAtomIndex(obj0,sele0);
      if(obj1)
        i1 = ObjectMoleculeGetAtomIndex(obj1,sele1);
      if(obj0&&obj1&&(i0>=0)&&(i1>=0)&&(obj0!=obj1)) {
        ObjectMoleculeVerifyChemistry(obj0,-1);
        ObjectMoleculeVerifyChemistry(obj1,-1);
        
        SelectorCreate(G,tmp_fuse_sele,NULL,obj0,1,NULL);
        sele2=SelectorIndexByName(G,tmp_fuse_sele);
        if(mode) {
          ObjectMoleculeOpRecInit(&op);
          op.code=OMOP_PrepareFromTemplate;
          op.ai=obj1->AtomInfo+i1;
          op.i1=mode;
          op.i2=0;
          op.i3=recolor;
          if(recolor)
            op.i4=obj1->Obj.Color;
          ExecutiveObjMolSeleOp(G,sele2,&op);
        }
        SelectorDelete(G,tmp_fuse_sele);

        if((obj0->AtomInfo[i0].protons==1)&&
           (obj1->AtomInfo[i1].protons==1))
          ObjectMoleculeFuse(obj1,i1,obj0,i0,0,move_flag);
        else if((obj0->AtomInfo[i0].protons!=1)&&
                (obj1->AtomInfo[i1].protons!=1))
          ObjectMoleculeFuse(obj1,i1,obj0,i0,1,move_flag);
        else 
          ErrMessage(G,"Fuse","Can't fuse between a hydrogen and a non-hydrogen");
      }
    }
  }
}

/*========================================================================*/
void ExecutiveSpheroid(PyMOLGlobals *G,char *name,int average)  /* EXPERIMENTAL */
{
  register CExecutive *I = G->Executive;
  CObject *os=NULL;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;

  if(strlen(name)) {
    os=ExecutiveFindObjectByName(G,name);
    if(!os)
      ErrMessage(G," Executive","object not found.");
    else if(os->type!=cObjectMolecule) {
      ErrMessage(G," Executive","bad object type.");
      os=NULL;
    }
  }
  
  if(os||(!strlen(name))) { /* sort one or all */
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule)
          if((!os)||(rec->obj==os)) {
            obj =(ObjectMolecule*)rec->obj;
            ObjectMoleculeCreateSpheroid(obj,average);  
            ObjectMoleculeInvalidate(obj,cRepAll,cRepInvRep,-1);
          }
    }
    SceneChanged(G);
  }
} 
/*========================================================================*/
void ExecutiveRebuildAll(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  PRINTFD(G,FB_Executive)
    " ExecutiveRebuildAll: entered.\n"
    ENDFD;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      switch(rec->obj->type) {
      case cObjectMolecule:
        if(SettingGetGlobal_b(G,cSetting_defer_builds_mode))
          ObjectMoleculeInvalidate((ObjectMolecule*)rec->obj,cRepAll,cRepInvPurge,-1);
        else
          ObjectMoleculeInvalidate((ObjectMolecule*)rec->obj,cRepAll,cRepInvRep,-1);           
        break;
      case cObjectMeasurement:
        ObjectDistInvalidateRep((ObjectDist*)rec->obj,cRepAll);
        break;
      case cObjectSurface:
      case cObjectMesh:
      case cObjectSlice:
      case cObjectAlignment:
      case cObjectCGO:
        if(rec->obj->fInvalidate) {
          rec->obj->fInvalidate((CObject*)rec->obj,cRepAll,cRepInvAll,-1);
        }
        break;
      }
    }
  }
  SeqChanged(G);
  SceneChanged(G);
}
/*========================================================================*/
void ExecutiveRebuildAllObjectDist(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      if(rec->obj->type==cObjectMeasurement) {
        ObjectDistInvalidateRep((ObjectDist*)rec->obj,cRepAll);
      }
    }
  }
  SceneInvalidate(G);
}
/*========================================================================*/
void ExecutiveUndo(PyMOLGlobals *G,int dir)
{
  register CExecutive *I = G->Executive;
  CObject *o;
  ObjectMolecule *obj=NULL,*compObj;
  SpecRec *rec = NULL;

  o = ExecutiveGetLastObjectEdited(G);
  PRINTFB(G,FB_Executive,FB_Debugging)
    " ExecutiveUndo: last object %p\n",(void*)o
    ENDFB(G);
  if(o)
    if(o->type==cObjectMolecule)
      obj = (ObjectMolecule*)o;
  /* make sure this is still a real object */
  if(obj) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject)
        if(rec->obj->type==cObjectMolecule) {
          compObj=(ObjectMolecule*)rec->obj;
          if(obj==compObj) {
            ObjectMoleculeUndo(obj,dir);
            break;
          }
        }
    }
  }
  
}
/*========================================================================*/
void ExecutiveSort(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  ObjectMolecule *obj;
  SpecRec *rec = NULL;
  ObjectMoleculeOpRec op;
  int sele;
#if 1
  if((!name)||(!name[0])) 
    name = cKeywordAll;
                            
  {
    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int changed = false;
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        switch(rec->type) {
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec,rec,next)) {
            if((rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
              obj =(ObjectMolecule*)rec->obj;
              ObjectMoleculeSort(obj);
              changed = true;
              sele=SelectorIndexByName(G,rec->name);
              if(sele>=0) {
                ObjectMoleculeOpRecInit(&op);
                op.code=OMOP_INVA;
                op.i1=cRepAll; 
                op.i2=cRepInvRep;
                ExecutiveObjMolSeleOp(G,sele,&op);
              }
            }
          }
          break;
        case cExecSelection:
          sele=SelectorIndexByName(G,rec->name);
          if(sele>=0) {
            op.code=OMOP_Sort;
            ExecutiveObjMolSeleOp(G,sele,&op);
            ObjectMoleculeOpRecInit(&op);
            op.code=OMOP_INVA;
            op.i1=cRepAll; 
            op.i2=cRepInvRep;
            ExecutiveObjMolSeleOp(G,sele,&op);
            ObjectMoleculeOpRecInit(&op);
          }
          break;
        case cExecObject:
          if(rec->obj->type==cObjectMolecule) {
            obj =(ObjectMolecule*)rec->obj;
            ObjectMoleculeSort(obj);
            changed = true;
            sele=SelectorIndexByName(G,rec->name);
            if(sele>=0) {
              ObjectMoleculeOpRecInit(&op);
              op.code=OMOP_INVA;
              op.i1=cRepAll; 
              op.i2=cRepInvRep;
              ExecutiveObjMolSeleOp(G,sele,&op);
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
    if(changed)
      SceneChanged(G);
  }
#else
  {
    CObject *os=NULL;
    int all_obj = false;

    if(strlen(name)) {
      os=ExecutiveFindObjectByName(G,name);
      if(!os) {
        if(!WordMatchExact(G,cKeywordAll,name,true))
          ErrMessage(G," Executive","object not found.");
        else
          all_obj=true;
      } else if(os->type!=cObjectMolecule)
        ErrMessage(G," Executive","bad object type.");
    } else {
      all_obj = true;
    }
  
    if(os||all_obj) { /* sort one or all */
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject)
          if(rec->obj->type==cObjectMolecule)
            if((rec->obj==os)||all_obj) {
              obj =(ObjectMolecule*)rec->obj;
              ObjectMoleculeSort(obj);
              sele=SelectorIndexByName(G,rec->obj->Name);
              if(sele>=0) {
                ObjectMoleculeOpRecInit(&op);
                op.code=OMOP_INVA;
                op.i1=cRepAll; 
                op.i2=cRepInvRep;
                ExecutiveObjMolSeleOp(G,sele,&op);
              }
            }
      }
      SceneChanged(G);
    }
  }
#endif

}
/*========================================================================*/
void ExecutiveRemoveAtoms(PyMOLGlobals *G,char *s1,int quiet)
{
  int sele;
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  ObjectMoleculeOpRec op;
  int flag = false;

  sele=SelectorIndexByName(G,s1);
  if(sele>=0)
    {
      while(ListIterate(I->Spec,rec,next))
        {
          if(rec->type==cExecObject)
            {
              if(rec->obj->type==cObjectMolecule)
                {
                  ObjectMoleculeOpRecInit(&op);
                  op.code = OMOP_Remove;
                  op.i1 = 0;
                  obj=(ObjectMolecule*)rec->obj;
                  ObjectMoleculeVerifyChemistry(obj,-1); /* remember chemistry for later */
                  ObjectMoleculeSeleOp(obj,sele,&op);
                  if(op.i1) {
                    if(!quiet) {
                      PRINTFD(G,FB_Editor)
                        " ExecutiveRemove-Debug: purging %i of %i atoms in %s\n",
                        op.i1,obj->NAtom,obj->Obj.Name
                        ENDFD;
                    }
                    ObjectMoleculePurge(obj);
                    if(!quiet) {
                      PRINTFB(G,FB_Editor,FB_Actions)
                        " Remove: eliminated %d atoms in model \"%s\".\n",
                        op.i1,obj->Obj.Name 
                        ENDFB(G);
                    }
                    flag=true;
                  }
                }
            }
        }
    }
  /*  if(!flag) {
      ErrMessage(G,"Remove","no atoms removed.");
      }*/
}
/*========================================================================*/
void ExecutiveAddHydrogens(PyMOLGlobals *G,char *s1,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_AddHydrogens; /* 4 passes completes the job */
    ExecutiveObjMolSeleOp(G,sele1,&op);    
  }
}
/*========================================================================*/
void ExecutiveFixHydrogens(PyMOLGlobals *G,char *s1,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_FixHydrogens; 
    ExecutiveObjMolSeleOp(G,sele1,&op);    
  }
}
/*========================================================================*/
void ExecutiveFlag(PyMOLGlobals *G,int flag,char *s1,int action,int quiet)
{
  int sele1;
  OrthoLineType buffer;
  ObjectMoleculeOpRec op;
  
  sele1 = SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    switch(action) {
    case 0: op.code = OMOP_Flag; break;
    case 1: op.code = OMOP_FlagSet; break;
    case 2: op.code = OMOP_FlagClear; break;
    default:
      op.code = OMOP_Flag;
      break;
    }
    op.i1 = (((unsigned int)1)<<flag);
    op.i2 = ((unsigned int)0xFFFFFFFF - (((unsigned int)1)<<flag));
    op.i3 = 0;
    op.i4 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op);    
    if(Feedback(G,FB_Executive,FB_Actions)) {
      if(!quiet) {
        switch(action) {
        case 0:
          if(op.i3) {
            PRINTF " Flag: flag %d is set in %d of %d atoms.\n", flag, op.i3, op.i4 ENDF(G);
          } else {
            PRINTF " Flag: flag %d cleared on all atoms.\n", flag ENDF(G);
          }
          break;
        case 1:
          PRINTF " Flag: flag %d set on %d atoms.\n", flag, op.i3 ENDF(G);
          break;
        case 2:
          PRINTF " Flag: flag %d cleared on %d atoms.\n", flag, op.i3 ENDF(G);
          break;
        }
      }
    }
    if((int)SettingGet(G,cSetting_auto_indicate_flags)) {
      sprintf(buffer,"(flag %d)",flag);
      SelectorCreate(G,cIndicateSele,buffer,NULL,true,NULL);
      ExecutiveSetObjVisib(G,cIndicateSele,true,false);
      SceneInvalidate(G);
    }
  }

}
/*========================================================================*/
float ExecutiveOverlap(PyMOLGlobals *G,char *s1,int state1,char *s2,int state2,float adjust)
{
  int sele1,sele2;
  float result=0.0;

  if(state1<0) state1=0;
  if(state2<0) state2=0;
                 
  sele1=SelectorIndexByName(G,s1);
  sele2=SelectorIndexByName(G,s2);

  if((sele1>=0)&&(sele2>=0))
    result = SelectorSumVDWOverlap(G,sele1,state1,sele2,state2,adjust);

  return(result);
}
/*========================================================================*/
void ExecutiveProtect(PyMOLGlobals *G,char *s1,int mode,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Protect;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op);    
    if(!quiet) {
      if(Feedback(G,FB_Executive,FB_Actions)) {
        if(op.i2) {
          if(mode) {
            PRINTF " Protect: %d atoms protected from movement.\n",op.i2 ENDF(G);
          } else {
            PRINTF " Protect: %d atoms deprotected.\n", op.i2 ENDF(G);
          }
        }
      }
    }
  }
}
/*========================================================================*/
void ExecutiveMask(PyMOLGlobals *G,char *s1,int mode,int quiet)
{
  int sele1;
  ObjectMoleculeOpRec op;
  
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_Mask;
    op.i1 = mode;
    op.i2 = 0;
    ExecutiveObjMolSeleOp(G,sele1,&op);
    if(!quiet) {
      if(Feedback(G,FB_Executive,FB_Actions)) {    
        if(op.i2) {
          if(mode) {
            PRINTF " Mask: %d atoms masked (cannot be picked or selected).\n",op.i2 ENDF(G);
          } else {
            PRINTF " Mask: %d atoms unmasked.\n", op.i2 ENDF(G);
          }
        }
      }
    }
    op.code = OMOP_INVA; /* need to invalidate all pickable representations */
    op.i1 = cRepAll;
    op.i2 = cRepInvPick;
    ExecutiveObjMolSeleOp(G,sele1,&op);    
  }
}
/*========================================================================*/
int ExecutiveStereo(PyMOLGlobals *G,int flag)
{
  int ok=1;
  int stereo_mode;

  switch(flag) {
  case -1:
    SettingSet(G,cSetting_stereo_shift,
               -SettingGet(G,cSetting_stereo_shift));
    /* shouldn't have to swap angle -- that's implicit
       SettingSet(cSetting_stereo_angle,-SettingGet(G,cSetting_stereo_angle));*/
    break;
  default:
    
    if(G->HaveGUI) {
      stereo_mode = (int)SettingGet(G,cSetting_stereo_mode);
      
      switch(stereo_mode) {
      case 1: /* hardware stereo-in-a-window*/
        SceneSetStereo(G,flag);
#ifndef _PYMOL_NOPY
        PSGIStereo(G,flag); /* does this have any effect anymore? */
#endif
        break;
      case 2: /* cross-eye stereo*/
      case 3: /* wall-eye */
      case 4: /* geo-wall */
      case 5: /* side-by-side */
        SceneSetStereo(G,flag);
        break;
      }
    }
  }
  SceneDirty(G);
  return(ok);
}
/*========================================================================*/
int ExecutiveRevalence(PyMOLGlobals *G,char *s1,char *s2,char *src,
                       int target_state,int source_state, int reset, int quiet)
{
  /*  register CExecutive *I=G->Executive; */
  int ok=true;
  int sele1,sele2;

  sele1=SelectorIndexByName(G,s1);
  sele2=SelectorIndexByName(G,s2);
  
  if((sele1>=0)&&(sele2>=0)) {
   if(src && (!src[0])) src=NULL;
   
   if(src) {
     int sele3=SelectorIndexByName(G,src);
     if(sele3>=0) {
       ObjectMolecule *obj3 = SelectorGetSingleObjectMolecule(G,sele3);
       if(!obj3) {
         ok=false;
         PRINTFB(G,FB_Editor,FB_Errors)
           "Editor-Warning: revalence can only source a single object at a time."
           ENDFB(G);
       } else {
         ObjectMoleculeOpRec op;
         
         ObjectMoleculeOpRecInit(&op);
         op.code = OMOP_RevalenceFromSource;
         op.i1 = sele1;
         op.i2 = sele2;
         op.i3 = target_state;
         op.obj3 = obj3;
         op.i4 = sele3;
         op.i5 = source_state;
         op.i6 = quiet;

         ExecutiveObjMolSeleOp(G,sele1,&op);
         
         /*
           if(ObjectMoleculeXferValences(obj1,sele1,sele2,target_state,obj3,sele3,source_state,quiet)) {
             ObjectMoleculeVerifyChemistry(obj1,target_state);
             ObjectMoleculeInvalidate(obj1,cRepAll,cRepInvBonds,target_state);
           }
         */

       }
     }
   } else { /* guess valences */
     ObjectMoleculeOpRec op;
     
     ObjectMoleculeOpRecInit(&op);
     op.code = OMOP_RevalenceByGuessing;
     op.i1 = sele1;
     op.i2 = sele2;
     op.i3 = target_state;
     op.i4 = reset;
     op.i6 = quiet;
     
     ExecutiveObjMolSeleOp(G,sele1,&op);
     
     }
  }
  return ok;
}

/*========================================================================*/
int ExecutiveBond(PyMOLGlobals *G,char *s1,char *s2,int order,int mode,int quiet)
{
  int ok=true;
  int sele1,sele2;
  int cnt;
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  int flag = false;

  sele1=SelectorIndexByName(G,s1);
  sele2=SelectorIndexByName(G,s2);
  
  if((sele1>=0)&&(sele2>=0)) {
    ObjectMolecule *obj1 = SelectorGetSingleObjectMolecule(G,sele1);
    ObjectMolecule *obj2 = SelectorGetSingleObjectMolecule(G,sele2);
    if((!obj1)||(!obj2)||(obj1!=obj2)) {
      if((!quiet)&&(mode==1)) {
        PRINTFB(G,FB_Editor,FB_Warnings)
          "Editor-Warning: bonds cannot be created between objects, only within.\n"
          ENDFB(G);
      }
    }
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          switch(mode) {
          case 1: /* add */
            cnt = ObjectMoleculeAddBond((ObjectMolecule*)rec->obj,sele1,sele2,order);
            if(cnt) {
              if(!quiet) {
                PRINTFB(G,FB_Editor,FB_Actions)
                  " Bond: %d bonds added to model \"%s\".\n",cnt,rec->obj->Name 
                  ENDFB(G);
                flag=true;
              }
            }
            break;
          case 2: /* adjust */
            cnt = ObjectMoleculeAdjustBonds((ObjectMolecule*)rec->obj,sele1,sele2,1,order);                    
            if(cnt) {
              if(!quiet) {
                PRINTFB(G,FB_Editor,FB_Actions)
                  " Valence: %d bond valences adjusted in model \"%s\".\n",cnt,rec->obj->Name 
                  ENDFB(G);
                flag=true;
              }
            }
            break;
          case 0: /* remove */
          default: 
            cnt = ObjectMoleculeRemoveBonds((ObjectMolecule*)rec->obj,sele1,sele2);
            if(cnt) {
              if(!quiet) {
                PRINTFB(G,FB_Editor,FB_Actions)
                  " Unbond: %d bonds removed from model \"%s\".\n",
                  cnt,rec->obj->Name 
                  ENDFB(G);
              }
              flag=true;
            }
          }
        }
      }
    }
    if(!flag) {
      if(!quiet) {
        switch(mode) {
        case 1:
          PRINTFB(G,FB_Editor,FB_Warnings) 
            "Bond-Warning: no bonds added."
            ENDFB(G);
          break;
        case 2:
          PRINTFB(G,FB_Editor,FB_Warnings) 
            "Valence-Warning: no bond valences changed."
            ENDFB(G);
          break;
        case 0:
        default:
          PRINTFB(G,FB_Editor,FB_Warnings) 
            "Unbond-Warning: no bonds removed."
            ENDFB(G);
          break;
        }
      }
    }
  } else if(sele1<0) {
    ok=ErrMessage(G,"ExecutiveBond","The first selection contains no atoms.");
  } else if(sele2<0) {
    ok=ErrMessage(G,"ExecutiveBond","The second selection contains no atoms.");
  }
  return ok;
}
/*========================================================================*/
int ExecutiveAngle(PyMOLGlobals *G,float *result, char *nam,
                   char *s1,char *s2, char *s3,int mode,
                   int labels,int reset,int zoom,int quiet,int state)
{
  int sele1,sele2,sele3;
  ObjectDist *obj;
  CObject *anyObj = NULL;
  sele1=SelectorIndexByName(G,s1);
  *result = 0.0F;
  if(!WordMatch(G,s2,cKeywordSame,true))
    sele2=SelectorIndexByName(G,s2);
  else {
    sele2 = sele1;
  }
  if(!WordMatch(G,s3,cKeywordSame,true))
    sele3=SelectorIndexByName(G,s3);
  else {
    sele3 = sele2;  
  }
  
  if((sele1>=0)&&(sele2>=0)&&(sele3>=0)) {
    anyObj = ExecutiveFindObjectByName(G,nam);
    if(anyObj) {
      if(anyObj->type!=cObjectMeasurement) {
        ExecutiveDelete(G,nam);
        anyObj=NULL;
      }
    }

    obj = ObjectDistNewFromAngleSele(G,(ObjectDist*)anyObj,
                                     sele1,sele2,sele3,
                                     mode,labels,result,reset,state);
    if(!obj) {
      if(!quiet)
        ErrMessage(G,"ExecutiveDistance","No angles found.");
    } else {
      *result = rad_to_deg(*result);
      if(!anyObj) {
        ObjectSetName((CObject*)obj,nam);
        ExecutiveManageObject(G,(CObject*)obj,zoom,quiet);
        ExecutiveSetRepVisib(G,nam,cRepLine,1);
        if(!labels)
          ExecutiveSetRepVisib(G,nam,cRepLabel,0);        
      }
    }
  } else if(sele1<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
               "The first selection contains no atoms.");
  } else if(sele2<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
               "The second selection contains no atoms.");
  } else if(sele3<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
               "The third selection contains no atoms.");
  }
  return(1);
}

/*========================================================================*/
int ExecutiveDihedral(PyMOLGlobals *G,float *result, char *nam,char *s1,
                        char *s2,char *s3,char *s4,int mode,
                        int labels,int reset,int zoom,int quiet,int state)
{
  int sele1,sele2,sele3,sele4;
  ObjectDist *obj;
  CObject *anyObj = NULL;
  sele1=SelectorIndexByName(G,s1);
  *result = 0.0F;
  if(!WordMatch(G,s2,cKeywordSame,true))
    sele2=SelectorIndexByName(G,s2);
  else {
    sele2 = sele1;
  }
  if(!WordMatch(G,s3,cKeywordSame,true))
    sele3=SelectorIndexByName(G,s3);
  else {
    sele3 = sele2;  
  }
  if(!WordMatch(G,s4,cKeywordSame,true))
    sele4=SelectorIndexByName(G,s4);
  else {
    sele4 = sele3;  
  }
  
  if((sele1>=0)&&(sele2>=0)&&(sele3>=0)&&(sele4>=0)) {
    anyObj = ExecutiveFindObjectByName(G,nam);
    if(anyObj) {
      if(anyObj->type!=cObjectMeasurement) {
        ExecutiveDelete(G,nam);
        anyObj=NULL;
      }
    }

    obj = ObjectDistNewFromDihedralSele(G,(ObjectDist*)anyObj,
                                        sele1,sele2,sele3,sele4,
                                        mode,labels,result,
                                        reset,state);
    if(!obj) {
      if(!quiet)
        ErrMessage(G,"ExecutiveDihedral","No angles found.");
    } else {
      *result = rad_to_deg(*result);
      if(!anyObj) {
        ObjectSetName((CObject*)obj,nam);
        ExecutiveManageObject(G,(CObject*)obj,zoom,quiet);
        ExecutiveSetRepVisib(G,nam,cRepLine,1);
        if(!labels)
          ExecutiveSetRepVisib(G,nam,cRepLabel,0);        
      }
    }
  } else if(sele1<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
                 "The first selection contains no atoms.");
  } else if(sele2<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
                 "The second selection contains no atoms.");
  } else if(sele3<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
                 "The third selection contains no atoms.");
  } else if(sele4<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
                 "The fourth selection contains no atoms.");
  }

  return 1;
}

int ExecutiveDist(PyMOLGlobals *G,float *result,char *nam,
                  char *s1,char *s2,int mode,float cutoff,
                  int labels,int quiet,int reset,
                  int state,int zoom)
{
  int sele1,sele2;
  ObjectDist *obj;
  CObject *anyObj = NULL;
  *result=0.0F;
  sele1=SelectorIndexByName(G,s1);
  if(!WordMatch(G,s2,"same",true))
    sele2=SelectorIndexByName(G,s2);
  else {
    sele2 = sele1;
  }
  if((sele1>=0)&&(sele2>=0)) {
    anyObj = ExecutiveFindObjectByName(G,nam);
    if(anyObj)
      if(reset || anyObj->type!=cObjectMeasurement) {
        ExecutiveDelete(G,nam);
        anyObj=NULL;
      }
    obj = ObjectDistNewFromSele(G,(ObjectDist*)anyObj,
                                sele1,sele2,mode,
                                cutoff,labels,
                                reset,result,state);
    if(!obj) {
      if(!quiet) 
        ErrMessage(G,"ExecutiveDistance",
                   "No such distances found.");
    } else {
      ObjectSetName((CObject*)obj,nam);
      ExecutiveManageObject(G,(CObject*)obj,zoom,quiet);
      ExecutiveSetRepVisib(G,nam,cRepLine,1);
      if(!labels)
        ExecutiveSetRepVisib(G,nam,cRepLabel,0);        
    }
  } else if(sele1<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
                 "The first selection contains no atoms.");
    if(reset)
      ExecutiveDelete(G,nam);
  } else if(sele2<0) {
    if(!quiet)
      ErrMessage(G,"ExecutiveDistance",
                 "The second selection contains no atoms.");
    if(reset)
      ExecutiveDelete(G,nam);
  }
  return 1;
}
/*========================================================================*/
float ExecutiveDistance(PyMOLGlobals *G,char *s1,char *s2)
{
  int sele1,sele2;
  float dist = -1.0;
  
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  
  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  sele1=SelectorIndexByName(G,s1);
  op1.i1=0;
  op2.i2=0;
  if(sele1>=0) {
    op1.code = OMOP_SUMC;
    op1.v1[0]=0.0;
    op1.v1[1]=0.0;
    op1.v1[2]=0.0;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
  } else {
    ErrMessage(G,"ExecutiveDistance","The first selection contains no atoms.");
  }
  
  sele2=SelectorIndexByName(G,s2);
  op2.i1=0;
  op2.i2=0;
  if(sele2>=0) {
    op2.code = OMOP_SUMC;
    op2.v1[0]=0.0;
    op2.v1[1]=0.0;
    op2.v1[2]=0.0;
    op2.i1=0;
    ExecutiveObjMolSeleOp(G,sele2,&op2);
  } else {
    ErrMessage(G,"ExecutiveDistance","The second selection contains no atoms.");
  }
  
  if(op1.i1&&op2.i1) {
    scale3f(op1.v1,1.0F/op1.i1,op1.v1);
    scale3f(op2.v1,1.0F/op2.i1,op2.v1);
    dist = (float)diff3f(op1.v1,op2.v1);
    PRINTFB(G,FB_Executive,FB_Results)
      " Distance: %8.3f [%i atom(s) to %i atom(s)]\n",
      dist,op1.i1,op2.i1
      ENDFB(G);
  } else {
    ErrMessage(G,"ExecutiveRMS","No atoms selected.");
  }
  return(dist);
}
/*========================================================================*/
char *ExecutiveNameToSeqAlignStrVLA(PyMOLGlobals *G,char *name,int state,int format,int quiet)
{
  char *result=NULL;    
  if((!name)||(!name[0])||(strcmp(name,"(all)")==0)) {
    /* use current alignment as the default */
    name = SettingGetGlobal_s(G,cSetting_seq_view_alignment);
    if(name[0]==0) {
      SpecRec *rec = NULL;
      register CExecutive *I = G->Executive;
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->visible) {
          if(rec->type==cExecObject)
            if(rec->obj->type==cObjectAlignment) {
              name = rec->obj->Name;              
              break;
            }
        }
      }
    }
  }
  if(!name) {
    ErrMessage(G," Executive","invalid alignment object name.");
  } else {
    CObject *obj=ExecutiveFindObjectByName(G,name);
    
    if(!obj) {
      ErrMessage(G," Executive","alignment object not found.");
    } else if(obj->type != cObjectAlignment) {
      ErrMessage(G," Executive","invalid object type.");
    } else {
      ObjectAlignmentAsStrVLA(G,(ObjectAlignment*)obj, state,format,&result);
    }
  }
  return(result);
}
/*========================================================================*/
char *ExecutiveSeleToPDBStr(PyMOLGlobals *G,char *s1,int state,int conectFlag,
                            int mode,char *ref_object,int ref_state,int quiet)
{
  char *result=NULL;
  ObjectMoleculeOpRec op1;
  int sele1;
  char end_str[] = "END\n";
  int model_count = 1;
  int actual_state = 0;
  int n_state = 1;
  int a;
  char model_record[50];
  int count=0,*counter=NULL;
  double matrix[16], inverse[16], *ref_mat = NULL;
  CObject *base = NULL;
  PDBInfoRec pdb_info;
  ObjectMolecule *obj = NULL;

  if(ref_object) {
    base=ExecutiveFindObjectByName(G,ref_object);
    if(base) {
      if(ref_state<-1) {
        ref_state = state;
      }
      if(ref_state<0) {
        ref_state = ObjectGetCurrentState(base,true);
      }
      if(ObjectGetTotalMatrix(base,ref_state,true,matrix)) {
        invert_special44d44d(matrix,inverse);
        ref_mat = inverse;
      }
    }
  }


  UtilZeroMem((void*)&pdb_info,sizeof(PDBInfoRec));
  ObjectMoleculeOpRecInit(&op1);
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    obj = SelectorGetSingleObjectMolecule(G,sele1);
    if(obj)
      if(obj->DiscreteFlag) {
        counter=&count; /* discrete objects need atom counters between states */
      }
  }
  op1.i2 = 0;
  op1.charVLA=VLAlloc(char,10000);

  if(state==-1) { /* multimodel PDB */
    n_state = ExecutiveCountStates(G,s1);
  }

  if(mode==1) {
    pdb_info.is_pqr_file = true;    
    pdb_info.pqr_workarounds = SettingGetGlobal_b(G,cSetting_pqr_workarounds);
  }

  for(a=0;a<n_state;a++) {
    switch(state) {
    case -1: /* multimodel */
      sprintf(model_record,"MODEL     %4d\n",model_count++);
      {
        ov_size len = op1.i2;
        UtilConcatVLA(&op1.charVLA,&len,model_record);
        op1.i2 = len;
      }
      actual_state = a;
      break;
    case -2: /* single state */
      actual_state=SceneGetState(G);

      if((actual_state!=0) && (sele1>=0) && SettingGetGlobal_b(G,cSetting_static_singletons)) {
        if(SelectorCountStates(G,sele1)==1) {
          actual_state = 0;
        }
      }
      break;
    default:
      actual_state = state;
      break;
    }
    
    if(conectFlag) {
      op1.i2=SelectorGetPDB(G,&op1.charVLA,op1.i2,sele1,
                            actual_state,conectFlag,&pdb_info,counter,ref_mat);
    } else {
      op1.i3 = 0; /* atIndex */
      if(sele1>=0) {
        op1.code = OMOP_PDB1;
        op1.i1 = actual_state;
        ExecutiveObjMolSeleOp(G,sele1,&op1);
      }
    }
    if((!(SettingGetGlobal_i(G,cSetting_pdb_no_end_record)))
       && !(pdb_info.is_pqr_file))
      /* terminate with END */
      {
        ov_size len = op1.i2;
        UtilConcatVLA(&op1.charVLA,&len,end_str);
        op1.i2 = len;
      }
    switch(state) {
    case -1:
      {
        ov_size len = op1.i2;
        UtilConcatVLA(&op1.charVLA,&len,"ENDMDL\n");
        op1.i2 = len;
      }
      break;
    }
  }

  /* terminate (just in case) */
  VLACheck(op1.charVLA,char,op1.i2+1);
  op1.charVLA[op1.i2]=0;
  op1.i2++;
  
  result=Alloc(char,op1.i2);
  memcpy(result,op1.charVLA,op1.i2);
  VLAFreeP(op1.charVLA);
  
  return(result);
}
/*========================================================================*/
PyObject *ExecutiveSeleToChemPyModel(PyMOLGlobals *G,char *s1,int state,
                                     char *ref_object,int ref_state)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;
  int sele1;
  CObject *base = NULL;
  double matrix[16], inverse[16], *ref_mat = NULL;

  if(ref_object) {
    base=ExecutiveFindObjectByName(G,ref_object);
    if(base) {
      if(ref_state<-1) {
        ref_state = state;
      }
      if(ref_state<0) {
        ref_state = ObjectGetCurrentState(base,true);
      }
      if(ObjectGetTotalMatrix(base,ref_state,true,matrix)) {
        invert_special44d44d(matrix,inverse);
        ref_mat = inverse;
      }
    }
  }

  sele1=SelectorIndexByName(G,s1);
  if(state<0) state=0;
  PBlock(G); /*   PBlockAndUnlockAPI();*/
  if(sele1>=0) {
    result=SelectorGetChemPyModel(G,sele1,state,ref_mat);
  }
  if(PyErr_Occurred()) PyErr_Print();
  PUnblock(G); /* PLockAPIAndUnblock();*/
  return(result);
#endif
}
/*========================================================================*/
int ExecutiveSeleToObject(PyMOLGlobals *G,char *name,char *s1,
                           int source,int target,
                          int discrete,int zoom, int quiet, int singletons)
{
  int sele1;
  int ok=false;
  ObjectNameType valid_name;

  UtilNCopy(valid_name, name, sizeof(valid_name));
  if(SettingGetGlobal_b(G,cSetting_validate_object_names)) {
    ObjectMakeValidName(valid_name);
    name = valid_name;
  }
  {
    int exists=(ExecutiveFindObjectMoleculeByName(G,name)!=NULL);
    
    sele1=SelectorIndexByName(G,s1);
    if(sele1>=0) {
      ok = SelectorCreateObjectMolecule(G,sele1,name,target,
                                        source,discrete,false,quiet,singletons);
      if(ok) {
        int sele2=SelectorIndexByName(G,name);
        ObjectMolecule *old_obj,*new_obj;
        old_obj = SelectorGetFirstObjectMolecule(G,sele1); /* get at least one object */
        new_obj = SelectorGetSingleObjectMolecule(G,sele2);
        
        /* first we need to make sure that the object being moved
           matches the target with respect to both the TTT and the
           object's state matrix (if any) */
        
        if(old_obj&&new_obj) {
          ExecutiveMatrixCopy(G,
                              old_obj->Obj.Name,
                              new_obj->Obj.Name, 
                              1, 1, /* TTT mode */
                              source,target,
                              false, 0, quiet);
          
          
          ExecutiveMatrixCopy(G,
                              old_obj->Obj.Name,
                              new_obj->Obj.Name, 
                              2, 2, /* Object state mode */
                              source,target,
                              false, 0, quiet);
          
          ExecutiveDoZoom(G,(CObject*)new_obj,!exists,zoom,true);        
        }
      }
    }
  }
  return ok;
}
/*========================================================================*/
void ExecutiveCopy(PyMOLGlobals *G,char *src,char *dst,int zoom)
{
  CObject *os;
  ObjectMolecule *oSrc,*oDst;
  SpecRec *rec1 = NULL,*rec2=NULL;
  int a;

  os=ExecutiveFindObjectByName(G,src);
  if(!os)
    ErrMessage(G," Executive","object not found.");
  else if(os->type!=cObjectMolecule)
    ErrMessage(G," Executive","bad object type.");
  else 
    {
      oSrc =(ObjectMolecule*)os;
      oDst = ObjectMoleculeCopy(oSrc);
      if(oDst) {
        strcpy(oDst->Obj.Name,dst);
        ExecutiveManageObject(G,(CObject*)oDst,zoom,false);
        rec1=ExecutiveFindSpec(G,oSrc->Obj.Name);
        rec2=ExecutiveFindSpec(G,oDst->Obj.Name);
        if(rec1&&rec2) {
          for(a=0;a<cRepCnt;a++)
            rec2->repOn[a]=rec1->repOn[a];
        }
        
        PRINTFB(G,FB_Executive,FB_Actions)
          " Executive: object %s created.\n",oDst->Obj.Name 
          ENDFB(G);
      }
    }
  SceneChanged(G);
}

/*========================================================================*/
void ExecutiveOrient(PyMOLGlobals *G,char *sele,double *mi,
                     int state,float animate,int complete,
                     float buffer,int quiet)
{
  double egval[3],egvali[3];
  double evect[3][3];
  float m[4][4],mt[4][4];
  float t[3];
  const float _0 = 0.0F;
  int a,b;

  if(!MatrixEigensolveC33d(G,mi,egval,egvali,(double*)evect)) {

	 normalize3d(evect[0]);
	 normalize3d(evect[1]);
	 normalize3d(evect[2]);

	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  m[a][b]=(float)evect[b][a]; /* fill columns */
		}
	 }

    for(a=0;a<3;a++) /* expand to 4x4 */
      {
        m[3][a]=0;
        m[a][3]=0;
      }
    m[3][3]=1.0;

    normalize3f(m[0]); /* cross normalization (probably unnec.)  */
    normalize3f(m[1]);
    normalize3f(m[2]);

    for(a=0;a<3;a++) /* convert to row-major */
      for(b=0;b<3;b++)
        mt[a][b]=m[b][a];

    cross_product3f(mt[0],mt[1],t);     /* insure right-handed matrix */
    if(dot_product3f(t,mt[2])<0.0) {
      mt[2][0] = -mt[2][0];
      mt[2][1] = -mt[2][1];
      mt[2][2] = -mt[2][2];
    }

    for(a=0;a<3;a++) /* convert back to column major */
      for(b=0;b<3;b++)
        m[a][b]=mt[b][a];

    if(animate<0.0F) {
      if(SettingGetGlobal_b(G,cSetting_animation))
        animate=SettingGetGlobal_f(G,cSetting_animation_duration);
      else
        animate=0.0F;
    }
    if(animate!=0.0F)
      ScenePrimeAnimation(G);

    {
      float old_mat[16];
      float new_mat[16];
      float x,y,z;
      copy44f(SceneGetMatrix(G),old_mat);
      
      SceneSetMatrix(G,m[0]); /* load matrix */
      
      /* there must  be a more elegant to get the PC on X and the SC
       * on Y then what is shown below, but I couldn't get it to work.
       * I tried swapping the eigen-columns around but either that is 
       * a bogus approach (?) or my code was buggy.  Hence the following...*/
      
      if((egval[0]<egval[2])&&(egval[2]<egval[1])) { /* X < Z < Y */
        SceneRotate(G,90,1,0,0); /*1<-->2*/
      } else if((egval[1]<egval[0])&&(egval[0]<egval[2])) { /* Y < X < Z */
        SceneRotate(G,90,0,0,1); /*0<-->1*/
      } else if((egval[1]<egval[2])&&(egval[2]<egval[0])) { /* Y < Z < X */
        SceneRotate(G,90,0,1,0); /*1<-->2*/
        SceneRotate(G,90,0,0,1); /*0<-->1*/
      } else if((egval[2]<egval[1])&&(egval[1]<egval[0])) { /* Z < Y < X */
        SceneRotate(G,90,0,1,0); /*0<-->2*/
      } else if((egval[2]<egval[0])&&(egval[0]<egval[1])) { /* Z < X < Y */
        SceneRotate(G,90,0,1,0); /*0<-->2*/
        SceneRotate(G,90,1,0,0); /*0<-->1*/
      }
      
      /* now choose orientation that has the least perturbation from the starting matrix */

      copy44f(SceneGetMatrix(G),new_mat);

      x = old_mat[0]*new_mat[0] + old_mat[4]*new_mat[4] + old_mat[ 8]*new_mat[ 8];
      y = old_mat[1]*new_mat[1] + old_mat[5]*new_mat[5] + old_mat[ 9]*new_mat[ 9];
      z = old_mat[2]*new_mat[2] + old_mat[6]*new_mat[6] + old_mat[10]*new_mat[10];
            
      if((x>_0)&&(y<_0)&&(z<_0)) {
        SceneRotate(G,180,1,0,0);
      } else if((x<_0)&&(y>_0)&&(z<_0)) {
        SceneRotate(G,180,0,1,0);        
      } else if((x<_0)&&(y<_0)&&(z>_0)) {
        SceneRotate(G,180,0,0,1);
      }
    }
    
    
    /* X < Y < Z  - do nothing - that's what we want */
    
    ExecutiveWindowZoom(G,sele,buffer,state,complete,false,quiet);
    if(animate!=0.0F)
      SceneLoadAnimation(G,animate,0);
    
  }
}
/*========================================================================*/
int ExecutiveLabel(PyMOLGlobals *G,char *s1,char *expr, int quiet, int eval)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  int cnt;

  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op1);
    op1.code = OMOP_LABL;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = eval;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    cnt = op1.i1;
    op1.code=OMOP_VISI;
    op1.i1=cRepLabel;
    op1.i2=1;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    op1.code = OMOP_INVA;
    op1.i1=cRepLabel; 
    op1.i2=cRepInvVisib;
    ExecutiveObjMolSeleOp(G,sele1,&op1);

    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Actions)
        " Label: labelled %i atoms.\n",cnt
        ENDFB(G);
    }
  } else {
    PRINTFB(G,FB_Executive,FB_Warnings)
      " Label: no atoms selected.\n"
      ENDFB(G);
  }
  return 1;
}
/*========================================================================*/
int ExecutiveIterate(PyMOLGlobals *G,char *s1,char *expr,int read_only,int quiet,PyObject *space)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRecInit(&op1);
  op1.i1=0;
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    op1.code = OMOP_ALTR;
    op1.s1 = expr;
    op1.i1 = 0;
    op1.i2 = read_only;
    op1.py_ob1 = space;
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    if(!quiet) {
      if(!read_only) {
        PRINTFB(G,FB_Executive,FB_Actions)
          " Alter: modified %i atoms.\n",op1.i1
          ENDFB(G);
      } else {
        PRINTFB(G,FB_Executive,FB_Actions)
          " Iterate: iterated over %i atoms.\n",op1.i1
          ENDFB(G);
      }
    }
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "ExecutiveIterate: No atoms selected.\n"
        ENDFB(G);
    }
  }
  return(op1.i1);
}
/*========================================================================*/
int ExecutiveSelectList(PyMOLGlobals *G,char *sele_name,char *s1,
                        int *list,int list_len,int state, int mode, int quiet)
{/* assumes a blocked Python interpreter */
  int ok=true;
  int n_eval=0;
  int sele0 = SelectorIndexByName(G,s1);
  int n_sele = 0;
  ObjectMolecule *obj = NULL;
  if(sele0>=0) obj = SelectorGetSingleObjectMolecule(G,sele0);
  if(obj) {
    int a;
    int index = 0;
    int check_state = true;
    CoordSet *cs = NULL;
    if(state==-2)
      state = SceneGetState(G);
    if(state==-3)
      state = ObjectGetCurrentState(&obj->Obj,true);
    if(state>=0) {
      cs = ObjectMoleculeGetCoordSet(obj,state);
    } else
      check_state = false;

    if(ok&&list) {
      if(list_len) {
        switch(mode) {
        case 0: /* object indices */
          for(a=0;a<list_len;a++) {
            list[a]--; /* convert 1-based indices to 0-based array offsets */
          }
          if(ok) 
            n_sele = SelectorCreateOrderedFromObjectIndices(G,sele_name,obj,list,list_len);
          break;
        case 1: /* atom identifier */
        case 2: /* rank */

          {
            OVOneToAny *o2a = OVOneToAny_New(G->Context->heap);
            AtomInfoType *ai;
            OVstatus res;
            OVreturn_word ret;
            int n_idx = 0;
            int *idx_list = VLAlloc(int,list_len);
            ai = obj->AtomInfo;

            for(a=0;a<obj->NAtom;a++) {
              ai->temp1 = -1;
              ai++;
            }

            /* create linked list using temp1 as "next" field */

            ai = obj->AtomInfo;
            for(a=0;a<obj->NAtom;a++) {
              if(mode==1) { /* id */
                index = ai[a].id;
              } else  /* rank */
                index = ai[a].rank;
              if((OVreturn_IS_ERROR( (res = OVOneToAny_SetKey(o2a,index, a))))) {
                if((OVreturn_IS_ERROR( (ret = OVOneToAny_GetKey(o2a,index))))) {
                  ok=false;
                } else {
                  int cur = ret.word;
                  while(1) {
                    if(ai[cur].temp1<0) {
                      ai[cur].temp1 = a; 
                      break; 
                    } else {
                      cur = ai[cur].temp1;
                    }
                  }
                }
              }
            }
            
            {
              int cur;
              for(a=0;a<list_len;a++) {
                index = list[a];
                if((OVreturn_IS_OK( (ret = OVOneToAny_GetKey(o2a,index))))) {
                  cur = ret.word;
                  while(cur>=0) {
                    if(check_state) {
                      if(cs) {
                        int ix;
                        if(obj->DiscreteFlag) {
                          if(cs==obj->DiscreteCSet[cur])
                            ix=obj->DiscreteAtmToIdx[a];
                          else
                            ix=-1;
                        } else 
                          ix=cs->AtmToIdx[a];
                        if(ix>=0) {
                          VLACheck(idx_list, int, n_idx);
                          idx_list[n_idx] = cur;
                          n_idx++;
                        }
                      }
                    } else {
                      VLACheck(idx_list, int, n_idx);
                      idx_list[n_idx] = cur;
                      n_idx++;
                    }
                    cur = ai[cur].temp1;
                  }
                }
              }
            }
            if(ok) 
              n_sele = SelectorCreateOrderedFromObjectIndices(G,sele_name,obj,idx_list,n_idx);
            OVOneToAny_DEL_AUTO_NULL(o2a);          
            VLAFreeP(idx_list);
          }
          break;
       
        }
      } else
        SelectorCreateEmpty(G,sele_name,true);
    }
  } else {
    PRINTFB(G,FB_Executive,FB_Errors)
      " SelectList-Error: selection cannot span more than one object.\n"
      ENDFB(G);
  }
  if(ok) {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Actions)
        " SelectList: modified %i atoms.\n",n_eval
        ENDFB(G);
    }
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "ExecutiveIterateList: An error occurred.\n"
        ENDFB(G);
    }
  }
  if(!ok)
    return -1;
  else
    return n_sele;
}


/*========================================================================*/
int ExecutiveIterateList(PyMOLGlobals *G,char *name,
                         PyObject *list,int read_only,int quiet,PyObject *space)
{
#ifdef _PYMOL_NOPY
  return -1;
#else
  int ok=true;
  int n_eval=0;
  int sele0 = SelectorIndexByName(G,name);
  PyObject *entry = NULL;
  ObjectMolecule *obj = NULL;
  if(sele0>=0) obj = SelectorGetSingleObjectMolecule(G,sele0);
  if(obj) {
    int n_atom = obj->NAtom;
    int list_len = 0;
    int a;
    int index = 0;
    char *expr = NULL;
    if(ok) ok=PyList_Check(list);
    if(ok) {
      list_len = PyList_Size(list);
      for(a=0;a<list_len;a++) {
        if(ok) entry=PyList_GetItem(list,a);
        if(ok) ok = PyList_Check(entry);
        if(ok) ok = (PyList_Size(entry)==2);
        if(ok) ok = PConvPyIntToInt(PyList_GetItem(entry,0),&index);
        if(ok) ok = PConvPyStrToStrPtr(PyList_GetItem(entry,1),&expr);
        if(ok) ok = ((index<=n_atom) && (index>0));
        if(ok) ok = PAlterAtom(G,obj->AtomInfo+index-1,expr,read_only,name,index-1,space);
        if(ok) n_eval++;
      }
    }
  } else {
    PRINTFB(G,FB_Executive,FB_Errors)
      " AlterList-Error: selection cannot span more than one object.\n"
      ENDFB(G);
  }
  if(ok) {
    if(!quiet) {
      if(!read_only) {
        PRINTFB(G,FB_Executive,FB_Actions)
          " AlterList: modified %i atoms.\n",n_eval
          ENDFB(G);
      } else {
        PRINTFB(G,FB_Executive,FB_Actions)
          " IterateList: iterated over %i atoms.\n",n_eval
          ENDFB(G);
      }
    }
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "ExecutiveIterateList: An error occurred.\n"
        ENDFB(G);
    }
  }
  if(!ok)
    return -1;
  else
    return n_eval;
#endif
}
/*========================================================================*/
void ExecutiveIterateState(PyMOLGlobals *G,int state,char *s1,char *expr,int read_only,
                           int atomic_props,int quiet,PyObject *space)
{
  int sele1;

  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    int start_state=0, stop_state=0;
    ObjectMoleculeOpRec op1;
    
    if(state>=0) {
      start_state = state;
      stop_state = state+1;
    } else {
      if((state==-2)||(state==-3)) { /* current state, TO DO: effective object state */
        state = SceneGetState(G);
        start_state = state;
        stop_state = state+1;
      } else if(state==-1) { /* all states (for the selection) */
        start_state = 0;
        stop_state = SelectorCountStates(G,sele1);
      }
    }
    ObjectMoleculeOpRecInit(&op1);
    op1.i1 = 0;

    for(state=start_state;state<stop_state;state++) {
      op1.code = OMOP_AlterState;
      op1.s1 = expr;
      op1.i2 = state;
      op1.i3 = read_only;
      op1.i4 = atomic_props;
      op1.py_ob1 = space;
      ExecutiveObjMolSeleOp(G,sele1,&op1);
    }
    if(!quiet) {
      if(!read_only) {
        PRINTFB(G,FB_Executive,FB_Actions)
          " AlterState: modified %i atom coordinate states.\n",op1.i1
          ENDFB(G);
      } else {
        PRINTFB(G,FB_Executive,FB_Actions)
          " IterateState: iterated over %i atom coordinate states.\n",op1.i1
          ENDFB(G);
      }
    }
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "ExecutiveIterateState: No atoms selected.\n"
        ENDFB(G);
    }
  }
}

typedef struct {
  int priority;
  float vertex[3];
  AtomInfoType *ai;
} FitVertexRec;

static int fVertexOrdered(FitVertexRec *array,int l, int r)
{
  return( array[l].priority <= array[r].priority );
}

static int fAtomOrdered(PyMOLGlobals *G,AtomInfoType **array,int l,int r)
{
  return(AtomInfoCompare(G, array[l], array[r]));
}

static int fAtomIDOrdered(AtomInfoType **array,int l,int r)
{
  return( array[l]->id <= array[r]->id );
}

static int fAtomRankOrdered(AtomInfoType **array,int l,int r)
{
  return( array[l]->rank <= array[r]->rank );
}

static int fAtomTemp1Ordered(AtomInfoType **array,int l,int r)
{
  return( array[l]->temp1 <= array[r]->temp1 );
}

static void PackSortedIndices(int n,int *x, int rec_size, void *data)
{
  register int a;
  for(a=0;a<n;a++) {
    if(a!=x[a]) {
    memcpy(((char*)data)+(a*rec_size),
           ((char*)data)+(x[a]*rec_size),
           rec_size);
    }
  }
}
  
/*========================================================================*/
int ExecutiveRMS(PyMOLGlobals *G,char *s1,char *s2,int mode,float refine,int max_cyc,
                   int quiet,char *oname,int state1,int state2,
                   int ordered_selections, int matchmaker, ExecutiveRMSInfo *rms_info)
{
  /* mode 0: measure rms without superposition
     mode 1: measure rms with trial superposition
     mode 2: measure rms with actual superposition */

  int sele1,sele2;
  float rms = -1.0;
  int a,b;
  float inv,*f,*f1,*f2;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType buffer;
  int *flag;
  int ok=true;
  int repeat;
  float v1[3],*v2;
  int matrix_mode = SettingGetGlobal_b(G,cSetting_matrix_mode);

  ObjectAlignment *align_to_update = NULL;

  sele1=SelectorIndexByName(G,s1);

  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  /* this function operates on stored coordinates -- thus transformation 
     matrices will need to be applied to the resulting atoms */

  if(sele1>=0) {
    if(state1<0) {
      op1.code = OMOP_AVRT;
    } else {
      op1.code = OMOP_StateVRT;
      op1.i1=state1;
    }
    op1.nvv1=0;
    op1.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1);
    if(mode==0) op1.i2 = true; /* if measuring current coordinates, then get global txfd values */
    if(matchmaker||(oname&&oname[0])) 
      op1.ai1VLA=(AtomInfoType**)VLAMalloc(1000,sizeof(AtomInfoType*),5,1);
    if(ordered_selections)
      op1.vp1=VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(G,sele1,&op1);
    for(a=0;a<op1.nvv1;a++)
      {
        inv=(float)op1.vc1[a]; /* average over coordinate sets */
        if(inv)
          {
            f=op1.vv1+(a*3);
            inv=1.0F/inv;
            *(f++)*=inv;
            *(f++)*=inv;
            *(f++)*=inv;
          }
      }
  }
  
  sele2=SelectorIndexByName(G,s2);
  if(sele2>=0) {

    if(state2<0) {
      op2.code = OMOP_AVRT;
    } else {
      op2.code = OMOP_StateVRT;
      op2.i1=state2;
    }
    op2.nvv1=0;
    op2.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
    op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1);
    if(mode==0) op2.i2 = true; /* if measuring current coordinates, then get global txfd values */
    if(matchmaker||(oname&&oname[0])) 
      op2.ai1VLA=(AtomInfoType**)VLAMalloc(1000,sizeof(AtomInfoType*),5,1);
    if(ordered_selections)
      op2.vp1=VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(G,sele2,&op2);
    for(a=0;a<op2.nvv1;a++)
      {
        inv=(float)op2.vc1[a]; /* average over coordinate sets */
        if(inv)
          {
            f=op2.vv1+(a*3);
            inv=1.0F/inv;
            *(f++)*=inv;
            *(f++)*=inv;
            *(f++)*=inv;
          }
      }
  }

  if(op1.vv1&&op2.vv1) {
    if(op1.nvv1&&op2.nvv1) {
      ObjectMolecule *mobile_obj = NULL;

      int n_pair = 0;

      if(! (mobile_obj = SelectorGetSingleObjectMolecule(G,sele1)) ) {
        if(mode!=2) {
          PRINTFB(G,FB_Executive,FB_Warnings)
            "Executive-Warning: Mobile selection spans more than one object.\n"
            ENDFB(G);
        } else {
          PRINTFB(G,FB_Executive,FB_Errors)
            "Executive-Error: Mobile selection spans more than one object. Aborting.\n"
            ENDFB(G);
          ok=false;
        }
      }
      
      if(ok && op1.nvv1 && op2.nvv1 && (matchmaker>0)) { /* matchmaker 0 is the default */
        int *idx1 = Alloc(int,op1.nvv1);
        int *idx2 = Alloc(int,op2.nvv1);
        int sort_flag = false;
        if(!(idx1&&idx2)) ok=false; else {
          switch(matchmaker) {
          case 1: /* by atom info */
            UtilSortIndexGlobals(G,op1.nvv1,op1.ai1VLA,idx1,(UtilOrderFnGlobals*)fAtomOrdered);
            UtilSortIndexGlobals(G,op2.nvv1,op2.ai1VLA,idx2,(UtilOrderFnGlobals*)fAtomOrdered);
            sort_flag = true;
            break;
          case 2: /* by matching atom identifiers */
            UtilSortIndex(op1.nvv1,op1.ai1VLA,idx1,(UtilOrderFn*)fAtomIDOrdered);
            UtilSortIndex(op2.nvv1,op2.ai1VLA,idx2,(UtilOrderFn*)fAtomIDOrdered);
            sort_flag = true;
            break;
          case 3: /* by matching atom ranks */
            UtilSortIndex(op1.nvv1,op1.ai1VLA,idx1,(UtilOrderFn*)fAtomRankOrdered);
            UtilSortIndex(op2.nvv1,op2.ai1VLA,idx2,(UtilOrderFn*)fAtomRankOrdered);
            sort_flag = true;
            break;
          case 4: /* by internal atom indexes (stored in temp1 kludge field) */
            UtilSortIndex(op1.nvv1,op1.ai1VLA,idx1,(UtilOrderFn*)fAtomTemp1Ordered);
            UtilSortIndex(op2.nvv1,op2.ai1VLA,idx2,(UtilOrderFn*)fAtomTemp1Ordered);
            sort_flag = true;
            break;
          }
          if(sort_flag) {  
            /* GOD this is SO ugly! */
            
            if(op1.vv1) {
              float *tmp = VLAlloc(float,op1.nvv1*3);
              if(!tmp) ok=false; else {
                UtilApplySortedIndices(op1.nvv1,idx1,3*sizeof(float),op1.vv1,tmp);
                VLAFreeP(op1.vv1);
                op1.vv1 = tmp;
              }
            }
            if(op1.vc1) {
              int *tmp = VLAlloc(int, op1.nvv1);
              if(!tmp) ok=false; else {
                UtilApplySortedIndices(op1.nvv1,idx1,sizeof(int),op1.vc1,tmp);
                VLAFreeP(op1.vc1);
                op1.vc1 = tmp;
              }
            }
            if(op1.vp1) {
              int *tmp = VLAlloc(int, op1.nvv1);
              if(!tmp) ok=false; else {
                UtilApplySortedIndices(op1.nvv1,idx1,sizeof(int),op1.vp1,tmp);
                VLAFreeP(op1.vp1);
                op1.vp1 = tmp;
              }
            }
            if(op1.ai1VLA) {
              AtomInfoType **tmp = VLAlloc(AtomInfoType*, op1.nvv1);
              if(!tmp) ok=false; else {
                UtilApplySortedIndices(op1.nvv1,idx1,sizeof(AtomInfoType*),op1.ai1VLA,tmp);
                VLAFreeP(op1.ai1VLA);
                op1.ai1VLA = tmp;
              }
            }

            if(op2.vv1) {
              float *tmp = VLAlloc(float,op2.nvv1*3);
              if(!tmp) ok=false; else {
                UtilApplySortedIndices(op2.nvv1,idx2,3*sizeof(float),op2.vv1,tmp);
                VLAFreeP(op2.vv1);
                op2.vv1 = tmp;
              }
            }
            if(op2.vc1) { 
              int *tmp = VLAlloc(int, op2.nvv1);
              if(!tmp) ok=false; else {
                UtilApplySortedIndices(op2.nvv1,idx2,sizeof(int),op2.vc1,tmp);
                VLAFreeP(op2.vc1);
                op2.vc1 = tmp;
              }
            }
            if(op2.vp1) {
              int *tmp = VLAlloc(int, op2.nvv1);
              if(!tmp) ok=false; else {
                UtilApplySortedIndices(op2.nvv1,idx2,sizeof(int),op2.vp1,tmp);
                VLAFreeP(op2.vp1);
                op2.vp1 = tmp;
              }
            }
            if(op2.ai1VLA) {
              AtomInfoType **tmp = VLAlloc(AtomInfoType*, op2.nvv1);
              if(!tmp) ok=false; else {
                UtilApplySortedIndices(op2.nvv1,idx2,sizeof(AtomInfoType*),op2.ai1VLA,tmp);
                VLAFreeP(op2.ai1VLA);
                op2.ai1VLA = tmp;
              }
            }

          }
        }

        if(matchmaker!=0) {
          int n1=0,n2=0,c1=0,c2=0;
          int cmp;
          
          while((n1<op1.nvv1)&&(n2<op2.nvv1)) {
            cmp = 0;
            switch(matchmaker) {
            case 1: /* insure that AtomInfoType matches */
              if(AtomInfoMatch(G,op1.ai1VLA[n1],op2.ai1VLA[n2]))
                cmp=0;
              else
                cmp=AtomInfoCompare(G,op1.ai1VLA[n1],op2.ai1VLA[n2]);
              printf("%d-%d %d-%d: %d\n",c1,n1,c2,n2,cmp);
              break;
            case 2: /* ID */
            case 3: /* rank */
              {
                int val1;
                int val2;
                
                switch(matchmaker) {
                case 2: /* ID */
                  val1=op1.ai1VLA[n1]->id;
                  val2=op2.ai1VLA[n2]->id;
                  break;
                case 3: /* rank */
                  val1=op1.ai1VLA[n1]->rank;
                  val2=op2.ai1VLA[n2]->rank;
                  break;
                case 4: /* index (via temp1) */
                  val1=op1.ai1VLA[n1]->temp1;
                  val2=op2.ai1VLA[n2]->temp1;
                  break;
                default:
                  val1=0;
                  val2=0;
                  break;
                }
                if(val1==val2)
                  cmp = 0;
                else if(val1<val2)
                  cmp = -1;
                else
                  cmp = 1;
              }
              break;
            }
            if(!cmp) { /* match found */
              idx1[c1++]=n1++;
              idx2[c2++]=n2++;
              n_pair++;
            } else if(cmp<0) { /* op1 below op2 */
              n1++;
            } else { /* op2 below op1 */
              n2++;
            }
          }

          if(n_pair) {
            if(op1.vv1) PackSortedIndices(n_pair,idx1,3*sizeof(float),op1.vv1);
            if(op1.vc1) PackSortedIndices(n_pair,idx1,sizeof(int),op1.vc1);
            if(op1.vp1) PackSortedIndices(n_pair,idx1,sizeof(int),op1.vp1);
            if(op1.ai1VLA) PackSortedIndices(n_pair,idx1,sizeof(AtomInfoType*),op1.ai1VLA);
            
            if(op2.vv1) PackSortedIndices(n_pair,idx2,3*sizeof(float),op2.vv1);
            if(op2.vc1) PackSortedIndices(n_pair,idx2,sizeof(int),op2.vc1);
            if(op2.vp1) PackSortedIndices(n_pair,idx2,sizeof(int),op2.vp1);
            if(op2.ai1VLA) PackSortedIndices(n_pair,idx2,sizeof(AtomInfoType*),op2.ai1VLA);
          }
        }
        FreeP(idx1);
        FreeP(idx2);
      } else if(op1.nvv1!=op2.nvv1) {
        sprintf(buffer,"Atom counts between selections don't match (%d vs %d)",
                op1.nvv1,op2.nvv1);
        ErrMessage(G,"ExecutiveRMS",buffer);
        n_pair = 0;
        ok=false;
      } else 
        n_pair = op1.nvv1;
      
      if(n_pair) {
        /* okay -- we're on track to do an alignment */

        if(ordered_selections&&op1.vp1&&op2.vp1) {
          /* if we expected ordered selections and have priorities, 
             then we may need to sort vertices */
          
          int sort_flag1 = false, sort_flag2 = false;
          int well_defined1 = true, well_defined2 = true;
          
          for(a=0;a<(n_pair-1);a++) {
            /*          printf("op1 vertex %d priority %d\n",a,op1.vp1[a]);
                        printf("op2 vertex %d priority %d\n",a,op2.vp1[a]);*/
            
            if(op1.vp1[a]>op1.vp1[a+1])
              sort_flag1 = true;
            else if(op1.vp1[a]==op1.vp1[a+1])
              well_defined1 = false;
            if(op2.vp1[a]>op2.vp1[a+1])
              sort_flag2 = true;
            else if(op2.vp1[a]==op2.vp1[a+1])
              well_defined2 = false;
          }
          
          if(sort_flag1||sort_flag2) {
            if(!(well_defined1||well_defined2)) {
              PRINTFB(G,FB_Executive,FB_Warnings) 
                "Executive-Warning: Ordering requested but not well defined.\n"
                ENDFB(G);
            } else {
              FitVertexRec *vert = Alloc(FitVertexRec,n_pair);

              if(sort_flag1) {
                float *src,*dst;
                src = op1.vv1;
                for(a=0;a<n_pair;a++) {              
                  vert[a].priority = op1.vp1[a];
                  dst=vert[a].vertex;
                  copy3f(src,dst);
                  src+=3;
                }
                UtilSortInPlace(G,vert,n_pair,sizeof(FitVertexRec),(UtilOrderFn*)fVertexOrdered);
                dst = op1.vv1;
                for(a=0;a<n_pair;a++) {              
                  src=vert[a].vertex;
                  copy3f(src,dst);
                  dst+=3;
                }
              }

              if(sort_flag2) {
                float *src,*dst;
                src = op2.vv1;
                for(a=0;a<n_pair;a++) {              
                  vert[a].priority = op2.vp1[a];
                  dst=vert[a].vertex;
                  copy3f(src,dst);
                  src+=3;
                }
                UtilSortInPlace(G,vert,n_pair,sizeof(FitVertexRec),(UtilOrderFn*)fVertexOrdered);
                dst = op2.vv1;
                for(a=0;a<n_pair;a++) {              
                  src=vert[a].vertex;
                  copy3f(src,dst);
                  dst+=3;
                }
              }
            
              FreeP(vert);
            }
          }
        }

        if(rms_info) {
          rms_info->initial_n_atom = n_pair;
          rms_info->n_cycles_run = 0;
          rms_info->final_n_atom = n_pair; /* in case there is no refinement */
        }
        

        if(mode!=0) {
          rms = MatrixFitRMSTTTf(G,n_pair,op1.vv1,op2.vv1,NULL,op2.ttt);
          if(rms_info) {
            rms_info->initial_rms = rms;
            rms_info->final_rms = rms;
          }
          repeat=true;
          b=0;
          while(repeat) {
            repeat=false;
            b++;
            if(b>max_cyc)
              break;
            if((refine>R_SMALL4)&&(rms>R_SMALL4)) {
              int n_next = n_pair;
              AtomInfoType **ai1,**ai2;
              
              flag=Alloc(int,n_pair);
              
              if(flag) {          
                for(a=0;a<n_pair;a++) {
                  MatrixTransformTTTfN3f(1,v1,op2.ttt,op1.vv1+(a*3));
                  v2=op2.vv1+(a*3);
                  if((diff3f(v1,v2)/rms)>refine) {
                    flag[a] = false;
                    repeat=true;
                  }
                  else
                    flag[a] = true;
                }
                f1 = op1.vv1;
                f2 = op2.vv1;
                ai1 = op1.ai1VLA;
                ai2 = op2.ai1VLA;
                for(a=0;a<n_pair;a++) {
                  if(!flag[a]) {
                    n_next--;
                  } else {
                    copy3f(op1.vv1+(3*a),f1);
                    copy3f(op2.vv1+(3*a),f2);
                    f1+=3;
                    f2+=3;
                    if(ai1&&ai2) { /* make sure we keep track of which atoms are aligned */
                      *(ai1++) = op1.ai1VLA[a];
                      *(ai2++) = op2.ai1VLA[a];
                    }
                  }
                }
                if(!quiet&&(n_next!=n_pair)) {
                  PRINTFB(G,FB_Executive,FB_Actions)
                    " ExecutiveRMS: %d atoms rejected during cycle %d (RMS=%0.2f).\n",n_pair-n_next,b,rms
                    ENDFB(G);
                }
                n_pair = n_next;
                FreeP(flag);
                if(n_pair) {
                  rms = MatrixFitRMSTTTf(G,n_pair,op1.vv1,op2.vv1,NULL,op2.ttt);            
                  if(rms_info) {
                    rms_info->n_cycles_run = b;
                    rms_info->final_n_atom = n_pair;
                    rms_info->final_rms = rms;
                  }
                } else
                  break;
              }
            }
          }
        } else { /* mode == 0 -- simple RMS, with no coordinate movement */
          rms = MatrixGetRMS(G,n_pair,op1.vv1,op2.vv1,NULL);
          if(rms_info) {
            rms_info->initial_rms = rms;
            rms_info->final_rms = rms;
          }
        }
      }
      if(!n_pair) {
        PRINTFB(G,FB_Executive,FB_Results) 
          " Executive: Error -- no atoms left after refinement!\n"
          ENDFB(G);
        ok = false;
      }

      if(ok) {
        if(!quiet) {
          PRINTFB(G,FB_Executive,FB_Results) 
            " Executive: RMS = %8.3f (%d to %d atoms)\n", rms,n_pair,n_pair
            ENDFB(G);
        }
        if(oname && oname[0]) {
#ifndef _PYMOL_1_x
          CGO *cgo = NULL;
          ObjectCGO *ocgo;
          int auto_save;
          
          cgo=CGONew(G);
          /*             CGOColor(cgo,1.0,1.0,0.0); 
                         CGOLinewidth(cgo,3.0);*/
          CGOBegin(cgo,GL_LINES);
          for(a=0;a<n_pair;a++) {
            CGOVertexv(cgo,op2.vv1+(a*3));
            MatrixTransformTTTfN3f(1,v1,op2.ttt,op1.vv1+(a*3));
            CGOVertexv(cgo,v1);
          }
          CGOEnd(cgo);
          CGOStop(cgo);
          ocgo = ObjectCGOFromCGO(G,NULL,cgo,0);
          ocgo->Obj.Color = ColorGetIndex(G,"yellow");
          ObjectSetName((CObject*)ocgo,oname);
          ExecutiveDelete(G,oname);
          auto_save = (int)SettingGet(G,cSetting_auto_zoom);
          SettingSet(G,cSetting_auto_zoom,0);
          ExecutiveManageObject(G,(CObject*)ocgo,-1,false);
          SettingSet(G,cSetting_auto_zoom,(float)auto_save);            
          SceneInvalidate(G);
#else
          {
            int align_state = state2;
            ObjectMolecule *trg_obj = SelectorGetSingleObjectMolecule(G,sele2);            

            if(align_state<0) {
              align_state = SceneGetState(G);
            }
            
            /* we're going to create/update an alignment object */
            
            {
              /* Get unique ids and construct the alignment vla */
              int *align_vla = VLAlloc(int, n_pair*3);
              
              {
                int *id_p = align_vla;
                int i;
                for(i=0;i<n_pair;i++) {
                  id_p[0] = AtomInfoCheckUniqueID(G,op2.ai1VLA[i]); /* target */
                  id_p[1] = AtomInfoCheckUniqueID(G,op1.ai1VLA[i]);
                  id_p[2] = 0;
                  id_p+=3;
                }
                VLASize(align_vla, int, n_pair*3);
              }
              {
                ObjectAlignment *obj = NULL;

                /* does object already exist? */
                {
                  CObject *execObj = ExecutiveFindObjectByName(G,oname);
                  if(execObj && (execObj->type != cObjectAlignment))
                    ExecutiveDelete(G,oname);
                  else
                    obj = (ObjectAlignment*)execObj;
                }
                obj = ObjectAlignmentDefine(G,obj,align_vla,align_state,true,trg_obj,mobile_obj);
                obj->Obj.Color = ColorGetIndex(G,"yellow");
                ObjectSetName((CObject*)obj,oname);
                ExecutiveManageObject(G,(CObject*)obj,0,quiet);
                align_to_update = obj;
                SceneInvalidate(G);
              }
              VLAFreeP(align_vla);
            }
          }
#endif
        }
        if(ok && mode==2) { 
          if(matrix_mode) {

            ObjectMolecule *src_obj,*trg_obj;
            src_obj = SelectorGetFirstObjectMolecule(G,sele1); /* get at least one object */
            trg_obj = SelectorGetSingleObjectMolecule(G,sele2);

            /* first we need to make sure that the object being moved
               matches the target with respect to both the TTT and the
               object's state matrix (if any) */

            if(src_obj&&trg_obj) {
              ExecutiveMatrixCopy(G,
                                  trg_obj->Obj.Name, 
                                  src_obj->Obj.Name,
                                  1, 1, /* TTT mode */
                                  state2,state1,
                                  false, 0, quiet);
              
              
              ExecutiveMatrixCopy(G,
                                  trg_obj->Obj.Name, 
                                  src_obj->Obj.Name,
                                  2, 2, /* Object state mode */
                                  state2,state1,
                                  false, 0, quiet);
              
              switch(matrix_mode) {
              case 1: /* TTTs */
                ExecutiveCombineObjectTTT(G,src_obj->Obj.Name,op2.ttt,true);
                break;
              case 2:
                {
                  double homo[16],*src_homo;
                  convertTTTfR44d(op2.ttt,homo);
                  if(ExecutiveGetObjectMatrix(G,src_obj->Obj.Name,state1,&src_homo,false)) {
                    left_multiply44d44d(src_homo,homo);
                    ExecutiveSetObjectMatrix(G,src_obj->Obj.Name,state1,homo);
                  }
                }
                break;
              }
              /* next we need to update the object's TTT matrix to reflect
                 the transformation */
            }
          } else { /* matrix_mode is zero -- legacy behavior */
            /* this will transform the actual coordinates */
            op2.code = OMOP_TTTF;
            ExecutiveObjMolSeleOp(G,sele1,&op2);
          }
        }
      }
    } else {
      ErrMessage(G,"ExecutiveRMS","No atoms selected.");
    }
  }

  if(align_to_update) {
    ObjectAlignmentUpdate(align_to_update);
  }
  
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  VLAFreeP(op1.vp1);
  VLAFreeP(op2.vp1);
  VLAFreeP(op1.ai1VLA);
  VLAFreeP(op2.ai1VLA);
  return(ok);
}
/*========================================================================*/
int *ExecutiveIdentify(PyMOLGlobals *G,char *s1,int mode)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  int *result = NULL;
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_Identify;
    op2.i1=0;
    op2.i1VLA=VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(G,sele1,&op2);
    result = op2.i1VLA;
    VLASize(result,int,op2.i1);
  } 
  return(result);
}
/*========================================================================*/
int ExecutiveIdentifyObjects(PyMOLGlobals *G,char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_IdentifyObjects;
    op2.obj1VLA=VLAlloc(ObjectMolecule*,1000);
    op2.i1VLA=VLAlloc(int,1000);
    op2.i1=0;
    ExecutiveObjMolSeleOp(G,sele1,&op2);
    VLASize(op2.i1VLA,int,op2.i1);
    VLASize(op2.obj1VLA,ObjectMolecule*,op2.i1);
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } 
  return(op2.i1);
}
/*========================================================================*/
int ExecutiveIndex(PyMOLGlobals *G,char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA)
{
  int sele1;
  ObjectMoleculeOpRec op2;
  sele1=SelectorIndexByName(G,s1);
  if(sele1>=0) {
    ObjectMoleculeOpRecInit(&op2);
    op2.code=OMOP_Index;
    op2.obj1VLA=VLAlloc(ObjectMolecule*,1000);
    op2.i1VLA=VLAlloc(int,1000);
    op2.i1=0;
    ExecutiveObjMolSeleOp(G,sele1,&op2);
    VLASize(op2.i1VLA,int,op2.i1);
    VLASize(op2.obj1VLA,ObjectMolecule*,op2.i1);
    (*indexVLA) = op2.i1VLA;
    (*objVLA) = op2.obj1VLA;
  } 
  return(op2.i1);
}
/*========================================================================*/
float *ExecutiveRMSStates(PyMOLGlobals *G,char *s1,int target,int mode,int quiet,int mix)
{
  int sele1;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  float *result = NULL;
  int ok=true;
  
  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.vv1=NULL;
  op2.vv1=NULL;
  sele1=SelectorIndexByName(G,s1);
  
  if(!SelectorGetSingleObjectMolecule(G,sele1)) {
    if(mode!=2) {
      PRINTFB(G,FB_Executive,FB_Warnings)
        "Executive-Warning: Mobile selection spans more than one object.\n"
        ENDFB(G);
    } else {
      PRINTFB(G,FB_Executive,FB_Errors)
        "Executive-Error: Mobile selection spans more than one object. Aborting.\n\n"
        ENDFB(G);
      ok=false;
    }
  }

  if(ok&&sele1>=0) {
    op1.code = OMOP_SVRT;
    op1.nvv1=0;
    op1.i1=target;
    op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    op1.i1VLA = VLAlloc(int,1000);
    ExecutiveObjMolSeleOp(G,sele1,&op1);

    op2.vv2=op1.vv1;
    op2.nvv2=op1.nvv1;
    op2.i1VLA=op1.i1VLA;
    op2.i2=target;
    op2.i1=mode;
    op2.i3 = mix;
    op2.f1VLA=VLAlloc(float,10);
    VLASize(op2.f1VLA,float,0); /* failsafe */
    op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,0);
    op2.code = OMOP_SFIT;
    op2.nvv1=0;
    ExecutiveObjMolSeleOp(G,sele1,&op2);
    result=op2.f1VLA;
    VLAFreeP(op1.vv1);
    VLAFreeP(op1.i1VLA);
    VLAFreeP(op2.vv1);
  } 
  return(result);
}
/*========================================================================*/
float ExecutiveRMSPairs(PyMOLGlobals *G,WordType *sele,int pairs,int mode)
{
  int sele1,sele2;
  int a,c;
  float rms=0.0,inv,*f;
  OrthoLineType buffer;

  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType combi,s1;

  ObjectMoleculeOpRecInit(&op1);
  ObjectMoleculeOpRecInit(&op2);
  op1.nvv1=0;
  op1.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
  op1.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1); /* auto-zero */
  op1.code = OMOP_AVRT;

  op2.nvv1=0;
  op2.vc1=(int*)VLAMalloc(1000,sizeof(int),5,1);
  op2.vv1=(float*)VLAMalloc(1000,sizeof(float),5,1); /* auto-zero */
  op2.code = OMOP_AVRT;

  strcpy(combi,"(");
  c=0;
  for(a=0;a<pairs;a++) {
    sele1=SelectorIndexByName(G,sele[c]);
    if(sele1>=0) ExecutiveObjMolSeleOp(G,sele1,&op1);
    strcat(combi,sele[c]);
    if(a<(pairs-1)) strcat(combi," or ");
    c++;
    sele2=SelectorIndexByName(G,sele[c]);
    if(sele2>=0) ExecutiveObjMolSeleOp(G,sele2,&op2);
    c++;
  }
  strcat(combi,")");
  for(a=0;a<op1.nvv1;a++)
    {
      inv=(float)op1.vc1[a];
      if(inv)
        {
          f=op1.vv1+(a*3);
          inv=1.0F/inv;
          *(f++)*=inv;
          *(f++)*=inv;
          *(f++)*=inv;
        }
    }
  for(a=0;a<op2.nvv1;a++)
    {
      inv=(float)op2.vc1[a];
      if(inv)
        {
          f=op2.vv1+(a*3);
          inv=1.0F/inv;
          *(f++)*=inv;
          *(f++)*=inv;
          *(f++)*=inv;
        }
    }
  if(op1.vv1&&op2.vv1) {
    if(op1.nvv1!=op2.nvv1) {
      sprintf(buffer,"Atom counts between selection sets don't match (%d != %d).",
              op1.nvv1,op2.nvv1);
      ErrMessage(G,"ExecutiveRMS",buffer);
    } else if(op1.nvv1) {
      if(mode!=0)
        rms = MatrixFitRMSTTTf(G,op1.nvv1,op1.vv1,op2.vv1,NULL,op2.ttt);
      else
        rms = MatrixGetRMS(G,op1.nvv1,op1.vv1,op2.vv1,NULL);
      PRINTFB(G,FB_Executive,FB_Results) 
        " ExecutiveRMS: RMS = %8.3f (%d to %d atoms)\n",
        rms,op1.nvv1,op2.nvv1
        ENDFB(G);

      op2.code = OMOP_TTTF;
      SelectorGetTmp(G,combi,s1);
      sele1=SelectorIndexByName(G,s1);
      ExecutiveObjMolSeleOp(G,sele1,&op2);
      SelectorFreeTmp(G,s1);
    } else {
      ErrMessage(G,"ExecutiveRMS","No atoms selected.");
    }
  }
  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  return(rms);
}
/*========================================================================*/
void ExecutiveUpdateObjectSelection(PyMOLGlobals *G,CObject *obj)
{
  if(obj->type==cObjectMolecule) {
    SelectorUpdateObjectSele(G,(ObjectMolecule*)obj);  
  }
}
/*========================================================================*/
int ExecutiveReset(PyMOLGlobals *G,int cmd,char *name)
{
  int ok=true;
  CObject *obj;
  if(!name[0]) {
    SceneResetMatrix(G);
    ExecutiveWindowZoom(G,cKeywordAll,0.0,-1,0,0,true); /* reset does all states */
  } else {
    obj = ExecutiveFindObjectByName(G,name);
    if(!obj)
      ok=false;
    else {
      ObjectResetTTT(obj);
      if(obj->fInvalidate)
        obj->fInvalidate(obj,cRepNone,cRepInvExtents,-1);

      SceneInvalidate(G);
    }
  }
  return(ok);
}
/*========================================================================*/
void ExecutiveDrawNow(PyMOLGlobals *G) 
{
  CExecutive *I = G->Executive;
  PRINTFD(G,FB_Executive)
    " ExecutiveDrawNow: entered.\n"
    ENDFD;

  if(PyMOL_GetIdleAndReady(G->PyMOL))
    OrthoExecDeferred(G);

  if(!SettingGet(G,cSetting_suspend_updates)) {

    if(G->HaveGUI && G->ValidContext) {
      glMatrixMode(GL_MODELVIEW); /* why is this necessary?  is it? */
    }

    ExecutiveUpdateSceneMembers(G);
    SceneUpdate(G,false);
    if(WizardUpdate(G))
      SceneUpdate(G,false);

    if(SettingGetGlobal_i(G,cSetting_stereo_mode)==4) {

      int width =  G->Option->winX;
      int height = G->Option->winY;
      
      glViewport(0,0,width/2,height);
      OrthoDoDraw(G,1);
      OrthoDoDraw(G,2);
      glViewport(0,0,width,height);

    } else {
      OrthoDoDraw(G,0);
    }

    if(G->HaveGUI && G->ValidContext) {
      if(I->CaptureFlag) {
        I->CaptureFlag = false;
        SceneCaptureWindow(G);
      }
    }
    PyMOL_NeedSwap(G->PyMOL);
  }

  PRINTFD(G,FB_Executive)
    " ExecutiveDrawNow: leaving.\n"
    ENDFD;
}
/*========================================================================*/
int ExecutiveCountStates(PyMOLGlobals *G,char *s1)
{
  register CExecutive *I = G->Executive;
  int sele1;
  int result = 0;
  int n_state;
  SpecRec *rec = NULL;
#if 1
  if((!s1)||(!s1[0])) s1 = cKeywordAll;
  {
    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G,s1,true,true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        switch(rec->type) {
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec,rec,next)) {
            if(rec->type==cExecObject) {
              if(rec->obj->fGetNFrame) {
                n_state = rec->obj->fGetNFrame(rec->obj);
                if(result<n_state)
                  result=n_state;
              }
            }
          }
          break;
        case cExecSelection:
          sele1=SelectorIndexByName(G,rec->name);
          if(sele1>=0) {
            SelectorUpdateTable(G,cSelectorUpdateTableAllStates,-1);
            n_state = SelectorGetSeleNCSet(G,sele1);
            if(result<n_state)
              result = n_state;
          }
          break;
        case cExecObject:
          if(rec->obj->fGetNFrame) {
            n_state = rec->obj->fGetNFrame(rec->obj);
            if(result<n_state)
              result=n_state;
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
#else

  if(s1)
    if(WordMatch(G,cKeywordAll,s1,true))
      s1 = NULL;
  if(!s1) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->fGetNFrame) {
          n_state = rec->obj->fGetNFrame(rec->obj);
          if(result<n_state)
            result=n_state;
        }
      }
    } 
  } else {
    sele1=SelectorIndexByName(G,s1);
    if(sele1>=0) {
      SelectorUpdateTable(G,cSelectorUpdateTableAllStates);
      result = SelectorGetSeleNCSet(G,sele1);
    }
  }
#endif
  return(result);

}
/*========================================================================*/
int ExecutiveRay(PyMOLGlobals *G,int width,int height,int mode,
                  float angle,float shift,int quiet,int defer, int antialias)
{
  if((mode==0) && 
     G->HaveGUI &&
     SettingGetGlobal_b(G,cSetting_auto_copy_images)) {
    /* force deferred behavior if copying image to clipboard */
    defer=1;
  }

  ExecutiveUpdateSceneMembers(G);

  if(defer && (mode==0)) {
    SceneDeferRay(G,width,height,mode,angle,shift,quiet,true,antialias);
  } else {
    SceneDoRay(G,width,height,mode,NULL,NULL,angle,shift,quiet,NULL,true,antialias);
  }
  return 1;
}
/*========================================================================*/
int *ExecutiveGetG3d(PyMOLGlobals *G)
{
  int *result = NULL;
  SceneRay(G,0,0,3,NULL,NULL,0.0F,0.0F,true,(G3dPrimitive**)&result,false,-1);
  return result;
}
/*========================================================================*/
int  ExecutiveSetBondSetting(PyMOLGlobals *G,int index,PyObject *tuple,
                             char *s1,char *s2,
                             int state,int quiet,int updates)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1,sele2;
  SettingName name;
  int unblock;
  int ok =true;
  int side_effects = false;
  int value_storage, *value_ptr;
  int value_type = 0;
  PRINTFD(G,FB_Executive)
    " ExecutiveSetBondSetting: entered. '%s' '%s'\n",s1,s2
    ENDFD;
  unblock = PAutoBlock(G);
  sele1 = SelectorIndexByName(G,s1);
  sele2 = SelectorIndexByName(G,s2);
  value_ptr = &value_storage;
  if((sele1>=0)&&(sele2>=0)) {
    int have_value=false;
    int type  = PyInt_AsLong(PyTuple_GetItem(tuple,0));
    PyObject *value = PyTuple_GetItem(tuple,1);
    if(value) {
      switch(type) {
      case cSetting_boolean:
        *(value_ptr) = PyInt_AsLong(PyTuple_GetItem(value,0));
        value_type = cSetting_boolean;
        have_value = true;
        break;
      case cSetting_int:
        *(value_ptr) = PyInt_AsLong(PyTuple_GetItem(value,0));
        value_type = cSetting_int;
        have_value = true;
        break;
      case cSetting_float:
        *(float*)value_ptr = (float)PyFloat_AsDouble(PyTuple_GetItem(value,0));
        value_type = cSetting_float;
        have_value = true;
        break;
      case cSetting_color:
        {
          int color_index=ColorGetIndex(G,PyString_AsString(PyTuple_GetItem(value,0)));
          if((color_index<0)&&(color_index>cColorExtCutoff))
            color_index = 0;
          *(value_ptr) = color_index;
          value_type = cSetting_color;
          have_value = true;
        }
        break;
      }
      if(have_value) {
        rec = NULL;
        while((ListIterate(I->Spec,rec,next))) {
          if((rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
            obj=(ObjectMolecule*)rec->obj;
            {
              int a,nBond = obj->NBond;
              int nSet = 0;
              BondType *bi = obj->Bond;
              AtomInfoType *ai1,*ai2, *ai = obj->AtomInfo;
              for(a=0;a<nBond;a++) {
                ai1 = ai + bi->index[0];
                ai2 = ai + bi->index[1];
                if((SelectorIsMember(G,ai1->selEntry,sele1) &&
                    SelectorIsMember(G,ai2->selEntry,sele2))||
                   (SelectorIsMember(G,ai2->selEntry,sele1) &&
                    SelectorIsMember(G,ai1->selEntry,sele2))) {
                  
                  int uid = AtomInfoCheckUniqueBondID(G,bi);
                  bi->has_setting = true;
                  SettingUniqueSetTypedValue(G,uid,index,value_type,value_ptr);
                  if(updates)
                    side_effects=true;
                  nSet++;
                }
                bi++;
              }
              if(nSet && !quiet) {
                SettingGetName(G,index,name);
                PRINTF
                  " Setting: %s set for %d bonds in object \"%s\".\n",
                  name,nSet, obj->Obj.Name
                  ENDF(G);
              }
            }
          }
        }
      }
    }
  }
  if(side_effects) {
    SettingGenerateSideEffects(G,index,s1,state);
    /*    SettingGenerateSideEffects(G,index,s2,state);*/
  }

  PAutoUnblock(G,unblock);
  return(ok);
#endif
}
/*========================================================================*/
int  ExecutiveUnsetBondSetting(PyMOLGlobals *G,int index,char *s1,char *s2,
                               int state,int quiet,int updates)
{
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  SettingName name;
  int unblock;
  int ok =true;
  int side_effects = false;
  int sele1,sele2;
  PRINTFD(G,FB_Executive)
    " ExecutiveSetSetting: entered. sele '%s' '%s'\n",s1,s2
    ENDFD;
  unblock = PAutoBlock(G);
  sele1 = SelectorIndexByName(G,s1);
  sele2 = SelectorIndexByName(G,s2);
  if((sele1>=0)&&(sele2>=0)) {
    rec = NULL;
    while((ListIterate(I->Spec,rec,next))) {
      if((rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
        obj=(ObjectMolecule*)rec->obj;
        {
          int a,nBond = obj->NBond;
          int nSet = 0;
          BondType *bi = obj->Bond;
          AtomInfoType *ai1,*ai2, *ai = obj->AtomInfo;
          for(a=0;a<nBond;a++) {
            ai1 = ai + bi->index[0];
            ai2 = ai + bi->index[1];
            if((SelectorIsMember(G,ai1->selEntry,sele1) &&
                SelectorIsMember(G,ai2->selEntry,sele2))||
               (SelectorIsMember(G,ai2->selEntry,sele1) &&
                SelectorIsMember(G,ai1->selEntry,sele2))) {
              int uid = AtomInfoCheckUniqueBondID(G,bi);
              bi->has_setting = true;
              SettingUniqueSetTypedValue(G,uid,index,cSetting_blank,NULL);
              if(updates)
                side_effects=true;
              nSet++;
            }
            bi++;
          }
          if(!quiet) {
            SettingGetName(G,index,name);
            PRINTF
              " Setting: %s unset for %d bonds in object \"%s\".\n",
              name, nSet, rec->obj->Name
              ENDF(G);
          }
        }
      }
    }
  }
  if(side_effects) {
    SettingGenerateSideEffects(G,index,s1,state);
    /*    SettingGenerateSideEffects(G,index,s2,state);*/
  }
  PAutoUnblock(G,unblock);
  return(ok);
}
/*========================================================================*/
int  ExecutiveSetSetting(PyMOLGlobals *G,int index,PyObject *tuple,char *sele,
                         int state,int quiet,int updates)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  OrthoLineType value;
  CSetting **handle=NULL;
  SettingName name;
  int nObj=0;
  int unblock;
  int ok =true;
  PRINTFD(G,FB_Executive)
    " ExecutiveSetSetting: entered. sele \"%s\"\n",sele
    ENDFD;
  unblock = PAutoBlock(G);
  if((!sele) || (sele[0]==0)) { /* global setting */
    ok = SettingSetFromTuple(G,NULL,index,tuple);
    if(ok) {
      if(!quiet) {
        if(Feedback(G,FB_Setting,FB_Actions)) {
          SettingGetTextValue(G,NULL,NULL,index,value);
          SettingGetName(G,index,name);
          PRINTF
            " Setting: %s set to %s.\n",name,value
            ENDF(G);
        }
      }
      if(updates) {
        SettingGenerateSideEffects(G,index,cKeywordAll,state);
      }
    }
  } 
#if 1
  else {
    int side_effects = false;

    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G,sele,true,true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        switch(rec->type) {
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec,rec,next)) {
            if(rec->type==cExecObject) {
              if(rec->obj->fGetSettingHandle) {
                handle = rec->obj->fGetSettingHandle(rec->obj,state);
                if(handle) {
                  SettingCheckHandle(G,handle);
                  ok = SettingSetFromTuple(G,*handle,index,tuple);
                  if(updates) 
                    side_effects = true;
                  nObj++;
                }
              }
            }
          }
          if(Feedback(G,FB_Setting,FB_Actions)) {
            if(nObj&&handle) {
              SettingGetTextValue(G,*handle,NULL,index,value);
              SettingGetName(G,index,name);
              if(!quiet) {
                if(state<0) {
                  PRINTF
                    " Setting: %s set to %s in %d objects.\n",name,value,nObj
                    ENDF(G);
                } else {
                  PRINTF
                    " Setting: %s set to %s in %d objects, state %d.\n",
                    name,value,nObj,state+1
                    ENDF(G);
                }
              }
            }
          }
          break;
        case cExecSelection:    
#if 1
          sele1 = SelectorIndexByName(G,rec->name);
          if(sele1>=0) {
            int have_atomic_value=false;
            int type  = PyInt_AsLong(PyTuple_GetItem(tuple,0));
            PyObject *value = PyTuple_GetItem(tuple,1);
            if(value) {
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_SetAtomicSetting;
              op.i1 = index;
              op.ii1 = &op.i3;
              switch(type) {
              case cSetting_boolean:
                *(op.ii1) = PyInt_AsLong(PyTuple_GetItem(value,0));
                op.i2 = cSetting_boolean;
                have_atomic_value = true;
                break;
              case cSetting_int:
                *(op.ii1) = PyInt_AsLong(PyTuple_GetItem(value,0));
                op.i2 = cSetting_int;
                have_atomic_value = true;
                break;
              case cSetting_float:
                *(float*)op.ii1 = (float)PyFloat_AsDouble(PyTuple_GetItem(value,0));
                op.i2 = cSetting_float;
                have_atomic_value = true;
                break;
              case cSetting_color:
                {
                  int color_index=ColorGetIndex(G,PyString_AsString(PyTuple_GetItem(value,0)));
                  if((color_index<0)&&(color_index>cColorExtCutoff))
                    color_index = 0;
                  *(op.ii1) = color_index;
                  op.i2 = cSetting_color;
                  have_atomic_value = true;
                }
                break;
              }
              if(have_atomic_value) {
                rec = NULL;
                while((ListIterate(I->Spec,rec,next))) {
                  if((rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
                    obj=(ObjectMolecule*)rec->obj;
                    op.i4 = 0;
                    ObjectMoleculeSeleOp(obj,sele1,&op);
                    if(op.i4) {
                      if(updates)
                        side_effects=true;
                      if(!quiet) {
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s set for %d atoms in object \"%s\".\n",
                            name,op.i4, rec->obj->Name
                          ENDF(G);
                      }
                    }
                  }
                }
              }
            }
          }
#else
          sele1 = SelectorIndexByName(G,rec->name);
          if(sele1>=0) {
            rec = NULL;
            while((ListIterate(I->Spec,rec,next))) {
              if((rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
                obj=(ObjectMolecule*)rec->obj;
                ObjectMoleculeOpRecInit(&op);
                op.code=OMOP_CountAtoms;
                op.i1=0;
                ObjectMoleculeSeleOp(obj,sele1,&op);
                if(op.i1&&rec->obj->fGetSettingHandle) {
                  handle = rec->obj->fGetSettingHandle(rec->obj,state);
                  if(handle) {
                    SettingCheckHandle(G,handle);
                    ok = SettingSetFromTuple(G,*handle,index,tuple);
                    if(ok) {
                      if(updates) 
                        side_effects = true;
                      if(!quiet) {
                        if(state<0) { /* object-specific */
                          if(Feedback(G,FB_Setting,FB_Actions)) {
                            SettingGetTextValue(G,*handle,NULL,index,value);
                            SettingGetName(G,index,name);
                            PRINTF
                              " Setting: %s set to %s in object \"%s\".\n",
                              name,value,rec->obj->Name
                              ENDF(G);
                          }
                        } else { /* state-specific */
                          if(Feedback(G,FB_Setting,FB_Actions)) {
                            SettingGetTextValue(G,*handle,NULL,index,value);
                            SettingGetName(G,index,name);
                            PRINTF
                              " Setting: %s set to %s in object \"%s\", state %d.\n",
                              name,value,rec->obj->Name,state+1
                              ENDF(G);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
#endif
          break;
        case cExecObject:
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
            if(handle) {
              SettingCheckHandle(G,handle);
              ok = SettingSetFromTuple(G,*handle,index,tuple);
              if(ok) {
                if(updates) 
                  side_effects = true;
                if(!quiet) {
                  if(state<0) { /* object-specific */
                    if(Feedback(G,FB_Setting,FB_Actions)) {
                      SettingGetTextValue(G,*handle,NULL,index,value);
                      SettingGetName(G,index,name);
                      PRINTF
                        " Setting: %s set to %s in object \"%s\".\n",
                        name,value,rec->obj->Name
                        ENDF(G);
                    }
                  } else { /* state-specific */
                    if(Feedback(G,FB_Setting,FB_Actions)) {
                      SettingGetTextValue(G,*handle,NULL,index,value);
                      SettingGetName(G,index,name);
                      PRINTF
                        " Setting: %s set to %s in object \"%s\", state %d.\n",
                        name,value,rec->obj->Name,state+1
                        ENDF(G);
                    }
                  }
                }
              }
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);

    if(side_effects)
      SettingGenerateSideEffects(G,index,sele,state);

  }
        
#else
  /* legacy code */

  else if(!strcmp(cKeywordAll,sele)) { /* all objects setting */
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject) {
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
            if(handle) {
              SettingCheckHandle(G,handle);
              ok = SettingSetFromTuple(G,*handle,index,tuple);
              nObj++;
            }
          }
        }
        if(nObj) {
          if(updates) 
            SettingGenerateSideEffects(G,index,sele,state);
        }
        if(Feedback(G,FB_Setting,FB_Actions)) {
          if(nObj&&handle) {
            SettingGetTextValue(G,*handle,NULL,index,value);
            SettingGetName(G,index,name);
            if(!quiet) {
              if(state<0) {
                PRINTF
                  " Setting: %s set to %s in %d objects.\n",name,value,nObj
                  ENDF(G);
              } else {
                PRINTF
                  " Setting: %s set to %s in %d objects, state %d.\n",
                  name,value,nObj,state+1
                  ENDF(G);
              }
            }
          }
        }
      }
  } else { /* based on a selection/object name */
    CObject *execObj = ExecutiveFindObjectByName(G,sele);
    if(execObj && 
       (execObj->type != cObjectMolecule))
      sele1 = -1;
    else
      sele1 = SelectorIndexByName(G,sele);
    while((ListIterate(I->Spec,rec,next)))
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule)
          {
            if(sele1>=0) {
              obj=(ObjectMolecule*)rec->obj;
              ObjectMoleculeOpRecInit(&op);
              op.code=OMOP_CountAtoms;
              op.i1=0;
              ObjectMoleculeSeleOp(obj,sele1,&op);
              if(op.i1&&rec->obj->fGetSettingHandle) {
                handle = rec->obj->fGetSettingHandle(rec->obj,state);
                if(handle) {
                  SettingCheckHandle(G,handle);
                  ok = SettingSetFromTuple(G,*handle,index,tuple);
                  if(ok) {
                    if(updates) 
                      SettingGenerateSideEffects(G,index,sele,state);
                    if(!quiet) {
                      if(state<0) { /* object-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetTextValue(G,*handle,NULL,index,value);
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s set to %s in object \"%s\".\n",
                            name,value,rec->obj->Name
                            ENDF(G);
                        }
                      } else { /* state-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetTextValue(G,*handle,NULL,index,value);
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s set to %s in object \"%s\", state %d.\n",
                            name,value,rec->obj->Name,state+1
                            ENDF(G);
                        }
                      }
                    }
                  }
                }
              }
            }
          } else if(strcmp(rec->obj->Name,sele)==0) {
            if(rec->obj->fGetSettingHandle) {
              handle = rec->obj->fGetSettingHandle(rec->obj,state);
              if(handle) {
                SettingCheckHandle(G,handle);
                ok = SettingSetFromTuple(G,*handle,index,tuple);
                if(ok) {
                  if(updates)
                    SettingGenerateSideEffects(G,index,sele,state);
                  if(!quiet) {
                    if(state<0) { /* object-specific */
                      if(Feedback(G,FB_Setting,FB_Actions)) {
                        SettingGetTextValue(G,*handle,NULL,index,value);
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s set to %s in object \"%s\".\n",
                          name,value,rec->obj->Name
                          ENDF(G);
                      }
                    } else { /* state-specific */
                      if(Feedback(G,FB_Setting,FB_Actions)) {
                        SettingGetTextValue(G,*handle,NULL,index,value);
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s set to %s in object \"%s\", state %d.\n",
                          name,value,rec->obj->Name,state+1
                          ENDF(G);
                      }
                    }
                  }
                }
              }
            }
          }
      }
  }
#endif

  PAutoUnblock(G,unblock);
  return(ok);
#endif
}

int  ExecutiveSetSettingFromString(PyMOLGlobals *G,
                                   int index,char *value,char *sele,
                                   int state,int quiet,int updates)
{
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  OrthoLineType value2;
  CSetting **handle=NULL;
  SettingName name;
  int nObj=0;
  int ok =true;

  PRINTFD(G,FB_Executive)
    " ExecutiveSetSetting: entered. sele \"%s\"\n",sele
    ENDFD;
  if(sele[0]==0) { /* global setting */
    ok = SettingSetFromString(G,NULL,index,value);
    if(ok) {
      if(!quiet) {
        if(Feedback(G,FB_Setting,FB_Actions)) {
          SettingGetTextValue(G,NULL,NULL,index,value2);
          SettingGetName(G,index,name);
          PRINTF
            " Setting: %s set to %s.\n",name,value2
            ENDF(G);
        }
      }
      if(updates) 
        SettingGenerateSideEffects(G,index,sele,state);
    }
  }
#if 1
  else {
    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G,sele,true,true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        switch(rec->type) {
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec,rec,next)) {
            if(rec->type==cExecObject) {
              if(rec->obj->fGetSettingHandle) {
                handle = rec->obj->fGetSettingHandle(rec->obj,state);
                if(handle) {
                  SettingCheckHandle(G,handle);
                  ok = SettingSetFromString(G,*handle,index,value);
                  if(updates) 
                    SettingGenerateSideEffects(G,index,rec->name,state);
                  nObj++;
                }
              }
            }
          }
          if(Feedback(G,FB_Setting,FB_Actions)) {
            if(nObj&&handle) {
              SettingGetTextValue(G,*handle,NULL,index,value2);
              SettingGetName(G,index,name);
              if(!quiet) {
                if(state<0) {
                  PRINTF
                    " Setting: %s set to %s in %d objects.\n",name,value2,nObj
                    ENDF(G);
                } else {
                  PRINTF
                    " Setting: %s set to %s in %d objects, state %d.\n",
                    name,value2,nObj,state+1
                    ENDF(G);
                }
              }
            }
          }
          break;
        case cExecSelection: 
#if 0
          /* this code has not yet been tested...*/

          sele1 = SelectorIndexByName(G,rec->name);
          if(sele1>=0) {
            int type;
            int value_store;
            if(SettingStringToTypedValue(G,index,value,&type,&value_store)) {
              ObjectMoleculeOpRecInit(&op);
              op.code = OMOP_SetAtomicSetting;
              op.i1 = index;
              op.i2 = type;
              op.ii1 = &value_store;
              rec = NULL;
              while((ListIterate(I->Spec,rec,next))) {
                if((rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
                  obj=(ObjectMolecule*)rec->obj;
                  op.i4 = 0;
                  ObjectMoleculeSeleOp(obj,sele1,&op);
                  if(op.i4) {
                    if(updates)
                      SettingGenerateSideEffects(G,index,rec->name,state);
                    if(!quiet) {
                      SettingGetName(G,index,name);
                      PRINTF
                        " Setting: %s set for %d atoms in object \"%s\".\n",
                        name,op.i4, rec->obj->Name
                        ENDF(G);
                    }
                  }
                }
              }
            }
          }
#else
          sele1 = SelectorIndexByName(G,rec->name);
          if(sele1>=0) {
            rec = NULL;
            while((ListIterate(I->Spec,rec,next))) {
              if( (rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
                obj=(ObjectMolecule*)rec->obj;
                ObjectMoleculeOpRecInit(&op);
                op.code=OMOP_CountAtoms;
                op.i1=0;
                ObjectMoleculeSeleOp(obj,sele1,&op);
                if(op.i1&&rec->obj->fGetSettingHandle) {
                  handle = rec->obj->fGetSettingHandle(rec->obj,state);
                  if(handle) {
                    SettingCheckHandle(G,handle);
                    ok = SettingSetFromString(G,*handle,index,value);
                    if(ok) {
                      if(updates) 
                        SettingGenerateSideEffects(G,index,sele,state);
                      if(!quiet) {
                        if(state<0) { /* object-specific */
                          if(Feedback(G,FB_Setting,FB_Actions)) {
                            SettingGetTextValue(G,*handle,NULL,index,value2);
                            SettingGetName(G,index,name);
                            PRINTF
                              " Setting: %s set to %s in object \"%s\".\n",
                              name,value2,rec->obj->Name
                              ENDF(G);
                          }
                        } else { /* state-specific */
                          if(Feedback(G,FB_Setting,FB_Actions)) {
                            SettingGetTextValue(G,*handle,NULL,index,value2);
                            SettingGetName(G,index,name);
                            PRINTF
                              " Setting: %s set to %s in object \"%s\", state %d.\n",
                              name,value2,rec->obj->Name,state+1
                              ENDF(G);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
#endif
          break;
        case cExecObject: 
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
            if(handle) {
              SettingCheckHandle(G,handle);
              ok = SettingSetFromString(G,*handle,index,value);
              if(ok) {
                if(updates)
                  SettingGenerateSideEffects(G,index,sele,state);
                if(!quiet) {
                  if(state<0) { /* object-specific */
                    if(Feedback(G,FB_Setting,FB_Actions)) {
                      SettingGetTextValue(G,*handle,NULL,index,value2);
                      SettingGetName(G,index,name);
                      PRINTF
                        " Setting: %s set to %s in object \"%s\".\n",
                        name,value2,rec->obj->Name
                        ENDF(G);
                    }
                  } else { /* state-specific */
                    if(Feedback(G,FB_Setting,FB_Actions)) {
                      SettingGetTextValue(G,*handle,NULL,index,value2);
                      SettingGetName(G,index,name);
                      PRINTF
                        " Setting: %s set to %s in object \"%s\", state %d.\n",
                        name,value2,rec->obj->Name,state+1
                        ENDF(G);
                    }
                  }
                }
              }
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
#else
  else if(!strcmp(cKeywordAll,sele)) { /* all objects setting */
    while(ListIterate(I->Spec,rec,next))
      {
        if(rec->type==cExecObject) {
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
            if(handle) {
              SettingCheckHandle(G,handle);
              ok = SettingSetFromString(G,*handle,index,value);
              nObj++;
            }
          }
        }
        if(nObj) {
          if(updates) 
            SettingGenerateSideEffects(G,index,sele,state);
        }
        if(Feedback(G,FB_Setting,FB_Actions)) {
          if(nObj&&handle) {
            SettingGetTextValue(G,*handle,NULL,index,value2);
            SettingGetName(G,index,name);
            if(!quiet) {
              if(state<0) {
                PRINTF
                  " Setting: %s set to %s in %d objects.\n",name,value2,nObj
                  ENDF(G);
              } else {
                PRINTF
                  " Setting: %s set to %s in %d objects, state %d.\n",
                  name,value2,nObj,state+1
                  ENDF(G);
              }
            }
          }
        }
      }
  } else { /* based on a selection/object name */
    CObject *execObj = ExecutiveFindObjectByName(G,sele);
    if(execObj && 
       (execObj->type != cObjectMolecule))
      sele1 = -1;
    else
      sele1 = SelectorIndexByName(G,sele);
    while((ListIterate(I->Spec,rec,next)))
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule)
          {
            if(sele1>=0) {
              obj=(ObjectMolecule*)rec->obj;
              ObjectMoleculeOpRecInit(&op);
              op.code=OMOP_CountAtoms;
              op.i1=0;
              ObjectMoleculeSeleOp(obj,sele1,&op);
              if(op.i1&&rec->obj->fGetSettingHandle) {
                handle = rec->obj->fGetSettingHandle(rec->obj,state);
                if(handle) {
                  SettingCheckHandle(G,handle);
                  ok = SettingSetFromString(G,*handle,index,value);
                  if(ok) {
                    if(updates) 
                      SettingGenerateSideEffects(G,index,sele,state);
                    if(!quiet) {
                      if(state<0) { /* object-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetTextValue(G,*handle,NULL,index,value2);
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s set to %s in object \"%s\".\n",
                            name,value2,rec->obj->Name
                            ENDF(G);
                        }
                      } else { /* state-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetTextValue(G,*handle,NULL,index,value2);
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s set to %s in object \"%s\", state %d.\n",
                            name,value2,rec->obj->Name,state+1
                            ENDF(G);
                        }
                      }
                    }
                  }
                }
              }
            }
          } else if(strcmp(rec->obj->Name,sele)==0) {
            if(rec->obj->fGetSettingHandle) {
              handle = rec->obj->fGetSettingHandle(rec->obj,state);
              if(handle) {
                SettingCheckHandle(G,handle);
                ok = SettingSetFromString(G,*handle,index,value);
                if(ok) {
                  if(updates)
                    SettingGenerateSideEffects(G,index,sele,state);
                  if(!quiet) {
                    if(state<0) { /* object-specific */
                      if(Feedback(G,FB_Setting,FB_Actions)) {
                        SettingGetTextValue(G,*handle,NULL,index,value2);
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s set to %s in object \"%s\".\n",
                          name,value2,rec->obj->Name
                          ENDF(G);
                      }
                    } else { /* state-specific */
                      if(Feedback(G,FB_Setting,FB_Actions)) {
                        SettingGetTextValue(G,*handle,NULL,index,value2);
                        SettingGetName(G,index,name);
                        PRINTF
                          " Setting: %s set to %s in object \"%s\", state %d.\n",
                          name,value2,rec->obj->Name,state+1
                          ENDF(G);
                      }
                    }
                  }
                }
              }
            }
          }
      }
  }
#endif
  return(ok);
}


int  ExecutiveSetObjSettingFromString(PyMOLGlobals *G,
                                      int index,char *value,CObject *obj,
                                      int state,int quiet,int updates)
{
  OrthoLineType value2;
  CSetting **handle=NULL;
  SettingName name;
  int ok =true;

  PRINTFD(G,FB_Executive)
    " ExecutiveSetObjSettingFromString: entered \n"
    ENDFD;
  if(!obj) { /* global */
    ok = SettingSetFromString(G,NULL,index,value);
    if(ok) {
      if(!quiet) {
        if(Feedback(G,FB_Setting,FB_Actions)) {
          SettingGetTextValue(G,NULL,NULL,index,value2);
          SettingGetName(G,index,name);
          PRINTF
            " Setting: %s set to %s.\n",name,value2
            ENDF(G);
        }
      }
      if(updates) 
        SettingGenerateSideEffects(G,index,obj->Name,state);
    }
  } 
  else { /* based on a single object */
    if(obj->fGetSettingHandle) {
      handle = obj->fGetSettingHandle(obj,state);
      if(handle) {
        SettingCheckHandle(G,handle);
        ok = SettingSetFromString(G,*handle,index,value);
        if(ok) {
          if(updates)
            SettingGenerateSideEffects(G,index,obj->Name,state);
          if(!quiet) {
            if(state<0) { /* object-specific */
              if(Feedback(G,FB_Setting,FB_Actions)) {
                SettingGetTextValue(G,*handle,NULL,index,value2);
                SettingGetName(G,index,name);
                PRINTF
                  " Setting: %s set to %s in object \"%s\".\n",
                  name,value2,obj->Name
                  ENDF(G);
              }
            } else { /* state-specific */
              if(Feedback(G,FB_Setting,FB_Actions)) {
                SettingGetTextValue(G,*handle,NULL,index,value2);
                SettingGetName(G,index,name);
                PRINTF
                  " Setting: %s set to %s in object \"%s\", state %d.\n",
                  name,value2,obj->Name,state+1
                  ENDF(G);
              }
            }
          }
        }
      }
    }
  }
  return(ok);
}
/*========================================================================*/
int  ExecutiveUnsetSetting(PyMOLGlobals *G,int index,char *sele,
                         int state,int quiet,int updates)
{
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;
  int sele1;
  ObjectMoleculeOpRec op;
  CSetting **handle=NULL;
  SettingName name;
  int nObj=0;
  int unblock;
  int ok =true;
  int side_effects = false;


  PRINTFD(G,FB_Executive)
    " ExecutiveSetSetting: entered. sele \"%s\"\n",sele
    ENDFD;
  unblock = PAutoBlock(G);
  if(sele[0]==0) { 
    /* do nothing -- in future, restore the default */
  } 
#if 1
  else {
    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G,sele,true,true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        switch(rec->type) {
        case cExecAll:
          rec = NULL;
          while(ListIterate(I->Spec,rec,next))
            {
              if(rec->type==cExecObject) {
                if(rec->obj->fGetSettingHandle) {
                  handle = rec->obj->fGetSettingHandle(rec->obj,state);
                  if(handle) {
                    SettingCheckHandle(G,handle);
                    ok = SettingUnset(*handle,index);
                    nObj++;
                  }
                }
              }
              if(nObj) {
                if(updates) 
                  side_effects = true;
              }
            }
          if(Feedback(G,FB_Setting,FB_Actions)) {
            if(nObj&&handle) {
              SettingGetName(G,index,name);
              if(!quiet) {
                if(state<0) {
                  PRINTF
                    " Setting: %s unset in %d objects.\n",name,nObj
                    ENDF(G);
                } else {
                  PRINTF
                    " Setting: %s unset in %d objects, state %d.\n",
                    name,nObj,state+1
                    ENDF(G);
                }
              }
            }
          }
          break;
        case cExecSelection:
#if 1
          sele1 = SelectorIndexByName(G,rec->name);
          if(sele1>=0) {
            ObjectMoleculeOpRecInit(&op);
            op.code = OMOP_SetAtomicSetting;
            op.i1 = index;
            op.i2 = cSetting_blank;
            op.ii1 = NULL;
            
            rec = NULL;
            while((ListIterate(I->Spec,rec,next))) {
              if((rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
                obj=(ObjectMolecule*)rec->obj;
                op.i4 = 0;
                ObjectMoleculeSeleOp(obj,sele1,&op);
                if(op.i4) {
                  if(updates)
                    side_effects=true;
                  if(!quiet) {
                    SettingGetName(G,index,name);
                    PRINTF
                      " Setting: %s unset for %d atoms in object \"%s\".\n",
                      name,op.i4, rec->obj->Name
                      ENDF(G);
                  }
                }
              }
            }
          }
#else
          sele1=SelectorIndexByName(G,rec->name);
          if(sele1>=0) {
            rec = NULL;
            while((ListIterate(I->Spec,rec,next))) {
              if((rec->type==cExecObject) && (rec->obj->type==cObjectMolecule)) {
                obj=(ObjectMolecule*)rec->obj;
                ObjectMoleculeOpRecInit(&op);
                op.code=OMOP_CountAtoms;
                op.i1=0;
                ObjectMoleculeSeleOp(obj,sele1,&op);
                if(op.i1&&rec->obj->fGetSettingHandle) {
                  handle = rec->obj->fGetSettingHandle(rec->obj,state);
                  if(handle) {
                    SettingCheckHandle(G,handle);
                    ok = SettingUnset(*handle,index);
                    if(ok) {
                      if(updates) 
                        side_effects = true;
                      if(!quiet) {
                        if(state<0) { /* object-specific */
                          if(Feedback(G,FB_Setting,FB_Actions)) {
                            SettingGetName(G,index,name);
                            PRINTF
                              " Setting: %s unset in object \"%s\".\n",
                              name,rec->obj->Name
                              ENDF(G);
                          }
                        } else { /* state-specific */
                          if(Feedback(G,FB_Setting,FB_Actions)) {
                            SettingGetName(G,index,name);
                            PRINTF
                              " Setting: %s unset in object \"%s\", state %d.\n",
                              name,rec->obj->Name,state+1
                              ENDF(G);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
#endif
          break;
        case cExecObject:
          if(rec->obj->fGetSettingHandle) {
            handle = rec->obj->fGetSettingHandle(rec->obj,state);
            if(handle) {
              SettingCheckHandle(G,handle);
              ok = SettingUnset(*handle,index);
              if(ok) {
                if(updates)
                  side_effects=true;
                if(!quiet) {
                  if(state<0) { /* object-specific */
                    if(Feedback(G,FB_Setting,FB_Actions)) {
                      SettingGetName(G,index,name);
                      PRINTF
                        " Setting: %s unset in object \"%s\".\n",
                        name,rec->obj->Name
                        ENDF(G);
                    }
                  } else { /* state-specific */
                    if(Feedback(G,FB_Setting,FB_Actions)) {
                      SettingGetName(G,index,name);
                      PRINTF
                        " Setting: %s unset in object \"%s\", state %d.\n",
                        name,rec->obj->Name,state+1
                        ENDF(G);
                    }
                  }
                }
              }
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
#else
  /* legacy */
  else {
    if(!strcmp(cKeywordAll,sele)) { /* all objects setting */
      
      while(ListIterate(I->Spec,rec,next))
        {
          if(rec->type==cExecObject) {
            if(rec->obj->fGetSettingHandle) {
              handle = rec->obj->fGetSettingHandle(rec->obj,state);
              if(handle) {
                SettingCheckHandle(G,handle);
                ok = SettingUnset(*handle,index);
                nObj++;
              }
            }
          }
          if(nObj) {
            if(updates) 
              side_effects = true;
          }
        }
      if(Feedback(G,FB_Setting,FB_Actions)) {
        if(nObj&&handle) {
          SettingGetName(G,index,name);
          if(!quiet) {
            if(state<0) {
              PRINTF
                " Setting: %s unset in %d objects.\n",name,nObj
                ENDF(G);
            } else {
              PRINTF
                " Setting: %s unset in %d objects, state %d.\n",
                name,nObj,state+1
                ENDF(G);
            }
          }
        }
      }
    } else { /* based on a selection/object name */
      sele1=SelectorIndexByName(G,sele);
      while((ListIterate(I->Spec,rec,next)))
        if(rec->type==cExecObject) {
          if(rec->obj->type==cObjectMolecule)
            {
              if(sele1>=0) {
                obj=(ObjectMolecule*)rec->obj;
                ObjectMoleculeOpRecInit(&op);
                op.code=OMOP_CountAtoms;
                op.i1=0;
                ObjectMoleculeSeleOp(obj,sele1,&op);
                if(op.i1&&rec->obj->fGetSettingHandle) {
                  handle = rec->obj->fGetSettingHandle(rec->obj,state);
                  if(handle) {
                    SettingCheckHandle(G,handle);
                    ok = SettingUnset(*handle,index);
                    if(ok) {
                      if(updates) 
                        side_effects = true;
                      if(!quiet) {
                        if(state<0) { /* object-specific */
                          if(Feedback(G,FB_Setting,FB_Actions)) {
                            SettingGetName(G,index,name);
                            PRINTF
                              " Setting: %s unset in object \"%s\".\n",
                              name,rec->obj->Name
                              ENDF(G);
                          }
                        } else { /* state-specific */
                          if(Feedback(G,FB_Setting,FB_Actions)) {
                            SettingGetName(G,index,name);
                            PRINTF
                              " Setting: %s unset in object \"%s\", state %d.\n",
                              name,rec->obj->Name,state+1
                              ENDF(G);
                          }
                        }
                      }
                    }
                  }
                }
              }
            } else if(strcmp(rec->obj->Name,sele)==0) {
              if(rec->obj->fGetSettingHandle) {
                handle = rec->obj->fGetSettingHandle(rec->obj,state);
                if(handle) {
                  SettingCheckHandle(G,handle);
                  ok = SettingUnset(*handle,index);
                  if(ok) {
                    if(updates)
                      side_effects=true;
                    if(!quiet) {
                      if(state<0) { /* object-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s unset in object \"%s\".\n",
                            name,rec->obj->Name
                            ENDF(G);
                        }
                      } else { /* state-specific */
                        if(Feedback(G,FB_Setting,FB_Actions)) {
                          SettingGetName(G,index,name);
                          PRINTF
                            " Setting: %s unset in object \"%s\", state %d.\n",
                            name,rec->obj->Name,state+1
                            ENDF(G);
                        }
                      }
                    }
                  }
                }
              }
            }
        }
    }
  }
#endif
  if(side_effects)
    SettingGenerateSideEffects(G,index,sele,state);
  PAutoUnblock(G,unblock);
  return(ok);
}
/*========================================================================*/
int ExecutiveColor(PyMOLGlobals *G,char *name,char *color,int flags,int quiet)
{
  /* flags: 
     0x1 -- ignore or suppress selection name matches
  */

  register CExecutive *I = G->Executive;
  int col_ind;
  int ok=false;
  col_ind = ColorGetIndex(G,color);
  if((!name) || (!name[0]))
    name = cKeywordAll;
  if(col_ind==-1) {
    ErrMessage(G,"Color","Unknown color.");
  } else {
    CTracker *I_Tracker= I->Tracker;
    SpecRec *rec = NULL;
    int n_atm=0;
    int n_obj=0;

#if 1
    int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        switch(rec->type) {
        case cExecSelection: 
        case cExecObject: 
        case cExecAll:
          if((rec->type==cExecSelection) || /* coloring a selection */
             (rec->type==cExecAll) || /* coloring all */
             ((rec->type==cExecObject) && /* coloring object and its backing selection */
              (rec->obj->type == cObjectMolecule))) {
            if(!(flags&0x1)) {
              int sele = SelectorIndexByName(G,rec->name);
              ObjectMoleculeOpRec op;
              if(sele>=0) {
                ok=true; 
                ObjectMoleculeOpRecInit(&op);
                op.code = OMOP_COLR;
                op.i1= col_ind;
                op.i2= n_atm;
                ExecutiveObjMolSeleOp(G,sele,&op);
                n_atm = op.i2;
                op.code=OMOP_INVA;
                op.i1=cRepAll; 
                op.i2=cRepInvColor;
                ExecutiveObjMolSeleOp(G,sele,&op);
              }
            }
          }
          break;
        }
        
        switch(rec->type) { /* sets object color */
        case cExecObject:
          rec->obj->Color=col_ind;
          if(rec->obj->fInvalidate)
            rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvColor,-1);
          n_obj++;
          ok=true;
          SceneInvalidate(G);
          break;
        case cExecAll:
          rec=NULL;
          while(ListIterate(I->Spec,rec,next)) {
            if(rec->type==cExecObject) {
              rec->obj->Color=col_ind;
              if(rec->obj->fInvalidate)
                rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvColor,-1);
              n_obj++;
              ok=true;
              SceneInvalidate(G);
            }
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);

#else
    /* per atom */
    if(!(flags&0x1)) {
      int sele = SelectorIndexByName(G,name);
      if(sele>=0) {
        ObjectMoleculeOpRec op;
        ok=true; 
        ObjectMoleculeOpRecInit(&op);
        op.code = OMOP_COLR;
        op.i1= col_ind;
        op.i2= 0;
        ExecutiveObjMolSeleOp(G,sele,&op);
        n_atm = op.i2;
        op.code=OMOP_INVA;
        op.i1=cRepAll; 
        op.i2=cRepInvColor;
        ExecutiveObjMolSeleOp(G,sele,&op);
      }
    }
    /* per object */
    if(strcmp(name,cKeywordAll)) {
      rec=ExecutiveFindSpec(G,name);
      if(rec) {
        if(rec->type==cExecObject) {
          rec->obj->Color=col_ind;
          if(rec->obj->fInvalidate)
            rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvColor,-1);
          n_obj++;
          ok=true;
          SceneInvalidate(G);
        }
      } 
    } else {
      rec=NULL;
      while(ListIterate(I->Spec,rec,next)) {
        if(rec->type==cExecObject) {
          rec->obj->Color=col_ind;
          if(rec->obj->fInvalidate)
            rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvColor,-1);
          n_obj++;
          ok=true;
          SceneInvalidate(G);
        }
      }
    }

#endif
    if(n_obj||n_atm) {
      char atms[]="s";
      char objs[]="s";
      if(n_obj<2) objs[0]=0;
      if(n_atm<2) atms[0]=0;
      if(!quiet) {

        if(n_obj&&n_atm) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Executive: Colored %d atom%s and %d object%s.\n",n_atm,atms,n_obj,objs
            ENDFB(G);
        } else if (n_obj) {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Executive: Colored %d object%s.\n",n_obj,objs
            ENDFB(G);
        } else {
          PRINTFB(G,FB_Executive,FB_Actions)
            " Executive: Colored %d atom%s.\n",n_atm,atms
            ENDFB(G);
        }
      }
    }
  }
  return(ok);
}
/*========================================================================*/
char *ExecutiveFindBestNameMatch(PyMOLGlobals *G,char *name)
{
  char *result;
  register CExecutive *I = G->Executive;
  SpecRec *rec=NULL,*best_rec = NULL;
  int best;
  int wm;

  best = 0;
  result = name;

  while(ListIterate(I->Spec,rec,next)) {
    wm = WordMatch(G,name,rec->name,true);
    if(wm<0) {
      best_rec = rec;
      best = wm;
      break;
    } else if ((best>0)&&(best<wm)) {
      best_rec=rec;
      best = wm;
    }
  }
  if(best_rec)
    result=best_rec->name;
  return(result);
}
/*========================================================================*/
static int count_objects(PyMOLGlobals *G,int public_only)
{
  int count = 0;
  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject) {
      if(!public_only)
        count++;
      else if(rec->obj->Name[0]!='_')
        count++;
    }
  }
  return count;
}

static SpecRec *ExecutiveFindSpec(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
#if 1
  { /* first, try for perfect, case-specific match */
    OVreturn_word result;
    if( OVreturn_IS_OK( (result = OVLexicon_BorrowFromCString(I->Lex,name)))) {
      if( OVreturn_IS_OK( (result = OVOneToOne_GetForward(I->Key, result.word)))) { 
        if(!TrackerGetCandRef(I->Tracker, result.word, (TrackerRef**)&rec)) {
          rec = NULL;
        }
      }
    }
    if(!rec) { /* otherwise try partial/case-nonspecific match */
      rec = ExecutiveAnyCaseNameMatch(G,name);
    }
  }
#else
  while(ListIterate(I->Spec,rec,next)) {
	 if(strcmp(rec->name,name)==0) 
		break;
  }
#endif
  return(rec);
}

/*========================================================================*/
void ExecutiveObjMolSeleOp(PyMOLGlobals *G,int sele,ObjectMoleculeOpRec *op) 
{
  register CExecutive *I=G->Executive;
  SpecRec *rec = NULL;
  ObjectMolecule *obj = NULL;

  if(sele>=0) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(rec->obj->type==cObjectMolecule) {
          obj=(ObjectMolecule*)rec->obj;
          ObjectMoleculeSeleOp(obj,sele,op);
        }
      }
    }
  }
}

/*========================================================================*/
int ExecutiveGetCameraExtent(PyMOLGlobals *G,char *name,float *mn,float *mx,int transformed,int state)
{
  int sele;
  ObjectMoleculeOpRec op;
  int flag = false;

  if((state==-2)||(state==-3)) /* TO DO: support per-object states */
    state=SceneGetState(G);

  PRINTFD(G,FB_Executive)
    " ExecutiveGetCameraExtent: name %s state %d\n",name,state
    ENDFD;
  
  sele=SelectorIndexByName(G,name);

  if(sele>=0) { 
    ObjectMoleculeOpRecInit(&op);
    if(state<0) {
      op.code = OMOP_CameraMinMax;
    } else {
      op.code = OMOP_CSetCameraMinMax;
      op.cs1 = state;
    }
	 op.v1[0]=FLT_MAX;
	 op.v1[1]=FLT_MAX;
	 op.v1[2]=FLT_MAX;
    op.v2[0]=-FLT_MAX;
    op.v2[1]=-FLT_MAX;
    op.v2[2]=-FLT_MAX;
    op.i1 = 0;
    op.i2 = transformed;
    op.mat1=SceneGetMatrix(G);

    ExecutiveObjMolSeleOp(G,sele,&op);

    PRINTFD(G,FB_Executive)
      " ExecutiveGetCameraExtent: minmax over %d vertices\n",op.i1
      ENDFD;
    if(op.i1)
      flag = true;
  }
  copy3f(op.v1,mn);
  copy3f(op.v2,mx);
  
  PRINTFD(G,FB_Executive)
    " ExecutiveGetCameraExtent: returning %d\n",flag
    ENDFD;

  return(flag);  
}

/*========================================================================*/
int ExecutiveGetExtent(PyMOLGlobals *G,char *name,float *mn,float *mx,
                       int transformed,int state,int weighted)
{
  int sele;
  ObjectMoleculeOpRec op,op2;
  register CExecutive *I=G->Executive;
  CObject *obj;
  int result = false;
  SpecRec *rec = NULL;
  float f1,f2,fmx;
  int a;

  if(WordMatch(G,cKeywordCenter,name,1)<0) {
    SceneGetPos(G,mn);
    copy3f(mn,mx);
    return 1;
  }
  if(WordMatch(G,cKeywordOrigin,name,1)<0) {
    SceneOriginGet(G,mn);
    copy3f(mn,mx);
    return 1;
  }

  PRINTFD(G,FB_Executive)
    " ExecutiveGetExtent: name %s state %d\n",name,state
    ENDFD;

  ObjectMoleculeOpRecInit(&op);
  ObjectMoleculeOpRecInit(&op2);  

  if((state==-2)||(state==-3)) { /* we want the currently displayed state */
    state=SceneGetState(G);
    op.include_static_singletons = true; /* make sure we get the static singletons too */
    op2.include_static_singletons = true;
  }

  op2.i1 = 0;
  op2.v1[0]=-1.0;
  op2.v1[1]=-1.0;
  op2.v1[2]=-1.0;
  op2.v2[0]=1.0;
  op2.v2[1]=1.0;
  op2.v2[2]=1.0;

  {
    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
    int have_atoms_flag = false;
    int have_extent_flag = false;

    /* first, compute atomic extents */

    if(weighted) {
      op2.i1 = 0;

      op2.v1[0]=0.0F;
      op2.v1[1]=0.0F;
      op2.v1[2]=0.0F;

      op.i1 = 0;

      op.v1[0]=FLT_MAX;
      op.v1[1]=FLT_MAX;
      op.v1[2]=FLT_MAX;

      op.v2[0]=-FLT_MAX;
      op.v2[1]=-FLT_MAX;
      op.v2[2]=-FLT_MAX;
    }

    /* first, handle molecular objects */

    {
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      
      while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
        if(rec) {
          switch(rec->type) {
          case cExecObject:
          case cExecSelection:
          case cExecAll:
            if(rec->type==cExecAll)
              sele=SelectorIndexByName(G,cKeywordAll);
            else
              sele=SelectorIndexByName(G,rec->name);            
            if(sele>=0) { 
              if(state<0) {
                op.code = OMOP_MNMX;
              } else {
                op.code = OMOP_CSetMinMax;
                op.cs1 = state;
              }
              op.i2 = transformed;
              ExecutiveObjMolSeleOp(G,sele,&op);
              if(op.i1) {
                have_atoms_flag = true;
              }
              PRINTFD(G,FB_Executive)
                " ExecutiveGetExtent: minmax over %d vertices\n",op.i1
                ENDFD;
            }
            
            if(weighted) {
              if(state<0) 
                op2.code = OMOP_SUMC;
              else {
                op2.code = OMOP_CSetSumVertices;
                op2.cs1 = state;
              }
              op2.i2 = transformed;
              ExecutiveObjMolSeleOp(G,sele,&op2);           
            }
            break;
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }
    if(have_atoms_flag)  have_extent_flag=true;

    /* now handle nonmolecular objects */

    {
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);

      while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
        if(rec) {
          switch(rec->type) {
          case cExecAll:
            rec = NULL;
            while(ListIterate(I->Spec,rec,next)) {
              if(rec->type==cExecObject) {
                obj=rec->obj;
                if(!obj->ExtentFlag) {
                  switch(obj->type) {
                  case cObjectMap:
                  case cObjectMesh:
                  case cObjectSurface:
                    if(!rec->obj->ExtentFlag) {
                      if(rec->obj->fUpdate) /* allow object to update extents, if necessary */
                        rec->obj->fUpdate(rec->obj);
                    }
                  }
                }
                if(obj->ExtentFlag) 
                  switch(obj->type) {
                  case cObjectMolecule:
                    break;
                    /* intentional fall-through */
                  default:
                    if(!have_extent_flag) {
                      copy3f(obj->ExtentMin,op.v1);
                      copy3f(obj->ExtentMax,op.v2);
                      have_extent_flag=true;
                    } else {
                      min3f(obj->ExtentMin,op.v1,op.v1);
                      max3f(obj->ExtentMax,op.v2,op.v2);
                    }
                    break;
                  }
              }
            }
            break;
          case cExecObject:
            obj=rec->obj;
            if(!obj->ExtentFlag) {
              switch(obj->type) {
              case cObjectMap:
              case cObjectMesh:
              case cObjectSurface:
                if(!rec->obj->ExtentFlag) {
                  if(rec->obj->fUpdate) /* allow object to update extents, if necessary */
                    rec->obj->fUpdate(rec->obj);
                }
              }
            }
            if(obj->ExtentFlag) 
              switch(obj->type) {
              case cObjectMolecule: /* will have been handled above... */
                break;
              default:
                if(!have_extent_flag) {
                  copy3f(obj->ExtentMin,op.v1);
                  copy3f(obj->ExtentMax,op.v2);
                  have_extent_flag=true;
                } else {
                  min3f(obj->ExtentMin,op.v1,op.v1);
                  max3f(obj->ExtentMax,op.v2,op.v2);
                }
                break;
              }
            break;
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }
   
    if(have_atoms_flag&&weighted) { 
      if(op2.i1) { 
        op2.v1[0]/=op2.i1; /* compute average */
        op2.v1[1]/=op2.i1;
        op2.v1[2]/=op2.i1;
        
        for (a=0;a<3;a++) { /* this puts origin at the weighted center */
          f1 = op2.v1[a] - op.v1[a];
          f2 = op.v2[a] - op2.v1[a];
          if(f1>f2) 
            fmx = f1;
          else
            fmx = f2;
          op.v1[a] = op2.v1[a] - fmx;
          op.v2[a] = op2.v1[a] + fmx;
        }
      }
    }

    if(have_extent_flag) {
      copy3f(op.v1,mn);
      copy3f(op.v2,mx);
    } else {
      zero3f(mn);
      zero3f(mx);
    }
    TrackerDelList(I_Tracker, list_id);
   
    result = have_extent_flag;

  }

  PRINTFD(G,FB_Executive)
    " ExecutiveGetExtent: returning %d\n",result
    ENDFD;

  return result;
}
/*========================================================================*/
static int ExecutiveGetMaxDistance(PyMOLGlobals *G,char *name,float *pos,float *dev,int transformed,int state)
{
  int sele;
  ObjectMoleculeOpRec op,op2;
  register CExecutive *I=G->Executive;
  CObject *obj;
  int flag = false;
  SpecRec *rec = NULL;
  float f1,fmx=0.0F;

  if((state==-2)||(state==-3)) /* TO DO: support per-object states */
    state=SceneGetState(G);

  PRINTFD(G,FB_Executive)
    " ExecutiveGetExtent: name %s state %d\n",name,state
    ENDFD;
  
  ObjectMoleculeOpRecInit(&op);
  ObjectMoleculeOpRecInit(&op2);

#if 1
  {
    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);

    op2.i1 = 0;
    op2.v1[0]=-1.0;
    op2.v1[1]=-1.0;
    op2.v1[2]=-1.0;
    op2.v2[0]=1.0;
    op2.v2[1]=1.0;
    op2.v2[2]=1.0;

    {
      /* first handle molecular objects */

      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      
      while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
        if(rec) {
          switch(rec->type) {
          case cExecObject:
          case cExecSelection:
          case cExecAll:
            if(rec->type==cExecAll)
              sele=SelectorIndexByName(G,cKeywordAll);
            else
              sele=SelectorIndexByName(G,rec->name);            
            if(sele>=0) { 
              if(state<0) {
                op.code = OMOP_MaxDistToPt;
              } else {
                op.code = OMOP_CSetMaxDistToPt;
                op.cs1 = state;
              }
              op.v1[0]=pos[0];
              op.v1[1]=pos[1];
              op.v1[2]=pos[2];
              op.i1 = 0;
              op.f1 = 0.0F;
              op.i2 = transformed;
              ExecutiveObjMolSeleOp(G,sele,&op);
              fmx = op.f1;
              
              if(op.i1)
                flag = true;
            }
            break;
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }

    {
      /* now handle nonmolecular objects */
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);

      while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
        if(rec) {
          switch(rec->type) {
          case cExecAll:
            rec = NULL;
            while(ListIterate(I->Spec,rec,next)) {
              if(rec->type==cExecObject) {
                obj=rec->obj;
                if(obj->ExtentFlag) {
                  switch(obj->type) {
                  case cObjectMolecule:
                    break;
                  default:
                    if(obj->ExtentFlag) {
                      f1 = (float)diff3f(obj->ExtentMin,pos);
                      if(fmx<f1) fmx = f1;
                      f1 = (float)diff3f(obj->ExtentMax,pos);
                      if(fmx<f1) fmx = f1;
                      flag = true;
                      break;
                    }
                  }
                }
              }
            }
            break;
          case cExecObject:
            obj=rec->obj;
            switch(rec->obj->type) {
            case cObjectMolecule:
              break;
            default:
              if(obj->ExtentFlag) {
                f1 = (float)diff3f(obj->ExtentMin,pos);
                if(fmx<f1) fmx = f1;
                f1 = (float)diff3f(obj->ExtentMax,pos);
                if(fmx<f1) fmx = f1;
                flag = true;
              }
              break;
            }
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);     
    }
    
    TrackerDelList(I_Tracker, list_id);
  }
#else
  {
  int all_flag = false;
  
    op2.i1 = 0;
    op2.v1[0]=-1.0;
    op2.v1[1]=-1.0;
    op2.v1[2]=-1.0;
    op2.v2[0]=1.0;
    op2.v2[1]=1.0;
    op2.v2[2]=1.0;
  
    if(WordMatch(G,cKeywordAll,name,true)<0) {
      all_flag=true;
    }
    sele=SelectorIndexByName(G,name);

    if(sele>=0) {
      if(state<0) {
        op.code = OMOP_MaxDistToPt;
      } else {
        op.code = OMOP_CSetMaxDistToPt;
        op.cs1 = state;
      }
      op.v1[0]=pos[0];
      op.v1[1]=pos[1];
      op.v1[2]=pos[2];
      op.i1 = 0;
      op.f1 = 0.0F;
      op.i2 = transformed;
      ExecutiveObjMolSeleOp(G,sele,&op);
      fmx = op.f1;

      if(op.i1)
        flag = true;
      if(all_flag) {
        while(ListIterate(I->Spec,rec,next)) {
          if(rec->type==cExecObject) {
            obj=rec->obj;
            if(obj->ExtentFlag) 
              switch(obj->type) {
              case cObjectMolecule:
                break;
              default:
                f1 = (float)diff3f(obj->ExtentMin,pos);
                if(fmx<f1) fmx = f1;
                f1 = (float)diff3f(obj->ExtentMax,pos);
                if(fmx<f1) fmx = f1;
                flag = true;
                break;
              }
          }
        }
      }
    } else {
      obj = ExecutiveFindObjectByName(G,name);
      if(obj) {
        switch(obj->type) {
        case cObjectMolecule:
          break;
        default:
          if(obj->ExtentFlag) {
            f1 = (float)diff3f(obj->ExtentMin,pos);
            if(fmx<f1) fmx = f1;
            f1 = (float)diff3f(obj->ExtentMax,pos);
            if(fmx<f1) fmx = f1;
            flag = true;
            break;
          }
        }
      } else if(all_flag) {
        rec=NULL;
        while(ListIterate(I->Spec,rec,next)) {
          if(rec->type==cExecObject) {
            obj=rec->obj;
            switch(rec->obj->type) {
            case cObjectMolecule:
              break;
            default:
              if(obj->ExtentFlag) {
                f1 = (float)diff3f(obj->ExtentMin,pos);
                if(fmx<f1) fmx = f1;
                f1 = (float)diff3f(obj->ExtentMax,pos);
                if(fmx<f1) fmx = f1;
              }
              break;
            }
          }
        }
      }
    }
  }
#endif
  *dev = fmx;
  return(flag);  
}
/*========================================================================*/
int ExecutiveWindowZoom(PyMOLGlobals *G,char *name,float buffer,
                        int state,int inclusive,float animate,
                        int quiet)
{
  float center[3],radius;
  float mn[3],mx[3],df[3];
  int sele0;
  int ok=true;


  PRINTFD(G,FB_Executive)
    " ExecutiveWindowZoom-DEBUG: entered\n"
    ENDFD;
  if(ExecutiveGetExtent(G,name,mn,mx,true,state,true)) {
    if(buffer!=0.0F) {
      buffer = buffer;
      mx[0]+=buffer;
      mx[1]+=buffer;
      mx[2]+=buffer;
      mn[0]-=buffer;
      mn[1]-=buffer;
      mn[2]-=buffer;
    }
    subtract3f(mx,mn,df);
    average3f(mn,mx,center);
    if(inclusive) {
      if(!ExecutiveGetMaxDistance(G,name,center,&radius,true,state))
        radius=0.0;
      radius+=buffer;
    } else {
      radius = df[0];
      if(radius<df[1]) radius=df[1];
      if(radius<df[2]) radius=df[2];
      radius=radius/2.0F;
    }
    if(radius<MAX_VDW) radius=MAX_VDW;
    PRINTFD(G,FB_Executive)
      " ExecutiveWindowZoom: zooming with radius %8.3f...state %d\n",radius,state
      ENDFD;
    PRINTFD(G,FB_Executive)
      " ExecutiveWindowZoom: on center %8.3f %8.3f %8.3f...\n",center[0],
      center[1],center[2]
      ENDFD;
    if(animate<0.0F) {
      if(SettingGetGlobal_b(G,cSetting_animation))
        animate=SettingGetGlobal_f(G,cSetting_animation_duration);
      else
        animate=0.0F;
    }
    if(animate!=0.0F)
      ScenePrimeAnimation(G);
    SceneOriginSet(G,center,false);
    SceneWindowSphere(G,center,radius);
    if(animate!=0.0F)
      SceneLoadAnimation(G,animate,0);
    else
      SceneAbortAnimation(G);
    SceneInvalidate(G);
  } else {

    sele0 = SelectorIndexByName(G,name);
    if(sele0>0) { /* any valid selection except "all" */
      /* no longer an error to zoom on an empty selection -- just has no effect */
      if(!quiet) {
        PRINTFB(G,FB_Executive, FB_Warnings)
          "ExecutiveWindowZoom-Warning: selection doesn't specify any coordinates.\n"
          ENDFB(G);
      }
    } else if(ExecutiveValidName(G,name)) {
      PRINTFD(G,FB_Executive)
        " ExecutiveWindowZoom-DEBUG: name valid, but no extents -- using default view\n"
        ENDFD;
      SceneSetDefaultView(G);
      SceneInvalidate(G);
    } else {
      ErrMessage(G,"ExecutiveWindowZoom","selection or object unknown.");
      ok=false;
    }
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveCenter(PyMOLGlobals *G,char *name,int state,
                    int origin,float animate, float *pos,int quiet)
{
  float center[3];
  float mn[3],mx[3],df[3];
  int sele0;
  int ok=true;
  int have_center = false;
  
  if(name && ExecutiveGetExtent(G,name,mn,mx,true,state,true)) {
    subtract3f(mx,mn,df);
    average3f(mn,mx,center);
    have_center = true;
    PRINTFD(G,FB_Executive)
      " ExecutiveCenter: centering state %d\n",state
      ENDFD;
    PRINTFD(G,FB_Executive)
      " ExecutiveCenter: on center %8.3f %8.3f %8.3f...\n",center[0],
      center[1],center[2]
      ENDFD;
  } else if(pos) {
    have_center = true;
    copy3f(pos,center);
  }
  if(have_center) {
    if(animate<0.0F) {
      if(SettingGetGlobal_b(G,cSetting_animation))
        animate=SettingGetGlobal_f(G,cSetting_animation_duration);
      else
        animate=0.0F;
    }

    if(animate!=0.0F)
      ScenePrimeAnimation(G);
    if(origin) 
      SceneOriginSet(G,center,false);
    SceneRelocate(G,center);
    SceneInvalidate(G);
    if(animate!=0.0F)
      SceneLoadAnimation(G,animate,0);
  } else {
    sele0 = SelectorIndexByName(G,name);
    if(sele0>=0) { /* any valid selection except "all" */
      if(!quiet) {
      /* no longer an error to center on an empty selection -- just have no effect */
        PRINTFB(G,FB_Executive, FB_Warnings)
          "ExecutiveCenter-Warning: selection doesn't specify any coordinates.\n"
          ENDFB(G);
      }
    } else if(ExecutiveValidName(G,name)) {
      SceneSetDefaultView(G);
      SceneInvalidate(G);
    } else {
      ErrMessage(G,"ExecutiveCenter","selection or object unknown.");
      ok=false;
    }
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveOrigin(PyMOLGlobals *G,char *name,int preserve,char *oname,float *pos,int state)
{
  float center[3];
  float mn[3],mx[3];
  int ok=true;
  CObject *obj = NULL;
  int have_center = false;
  if(oname && oname[0]) {
    obj = ExecutiveFindObjectByName(G,oname);
    if(!obj)
      ok=false;
  }
  if(ok) {
    if(name && name[0]) {
      ok = ExecutiveGetExtent(G,name,mn,mx,true,state,true);
      if(ok) {
        average3f(mn,mx,center);
        have_center = true;
      }
    } else if(pos) {
      copy3f(pos,center);
      have_center = true;
    }
  }
  if(ok && have_center) {
    if(obj) {
      ObjectSetTTTOrigin(obj,center);
      PRINTFB(G,FB_Executive,FB_Blather)
        " ExecutiveCenter: origin for %s set to %8.3f %8.3f %8.3f\n",
        oname,center[0],center[1],center[2]
        ENDFB(G);
    } else {
      PRINTFB(G,FB_Executive,FB_Blather)
        " ExecutiveCenter: scene origin set to %8.3f %8.3f %8.3f\n",
        center[0],center[1],center[2]
        ENDFB(G);
      SceneOriginSet(G,center,preserve);
    }
    SceneInvalidate(G);
  } else
    ok=false;
  return(ok);
}
/*========================================================================*/
int ExecutiveGetMoment(PyMOLGlobals *G,char *name,double *mi,int state)
{
  int sele;
  ObjectMoleculeOpRec op;
  int a,b;
  int c=0;

  if((state==-2)||(state==-3)) /* TO DO: support per-object states */
    state=SceneGetState(G);

  sele=SelectorIndexByName(G,name);
  if(sele>=0) {
    ObjectMoleculeOpRecInit(&op);
    if(state<0) {
      op.code = OMOP_SUMC;
    } else {
      op.code = OMOP_CSetSumVertices;
      op.cs1=state;
    }
    
    op.v1[0]=0.0;
    op.v1[1]=0.0;
    op.v1[2]=0.0;
    op.i1=0;
    op.i2=0; /* untransformed...is this right? */
	 
	 ExecutiveObjMolSeleOp(G,sele,&op);
	 
	 if(op.i1) { /* any vertices? */
		c+=op.i1;
		scale3f(op.v1,1.0F/op.i1,op.v1); /* compute raw average */
      if(state<0) {
        op.code = OMOP_MOME;		
      } else {
        op.code = OMOP_CSetMoment;
        op.cs1=state;
      }
		for(a=0;a<3;a++)
		  for(b=0;b<3;b++)
			 op.d[a][b]=0.0;
		ExecutiveObjMolSeleOp(G,sele,&op);			 
      {
        double *p = mi;
        for(a=0;a<3;a++)
          for(b=0;b<3;b++)
            *(p++)=op.d[a][b];
      }
	 }
  } else {
    identity33d(mi);
  }

  return(c);
}
/*========================================================================*/
int ExecutiveSetObjVisib(PyMOLGlobals *G,char *name,int onoff,int parents)
{
  register CExecutive *I = G->Executive;
  PRINTFD(G,FB_Executive)
    " ExecutiveSetObjVisib: entered.\n"
    ENDFD;
#if 1
  {
    CTracker *I_Tracker= I->Tracker;
    SpecRec *rec;
    int list_id = ExecutiveGetNamesListFromPattern(G,name,true,false);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    int suppress_hidden = SettingGetGlobal_b(G,cSetting_suppress_hidden);
    int hide_underscore = SettingGetGlobal_b(G,cSetting_hide_underscore_names);
    if(suppress_hidden && hide_underscore)
      ExecutiveUpdateGroups(G,false);
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {

      if(rec) {
        switch(rec->type) {
        case cExecAll:
          {
            SpecRec *tRec=NULL;
            while(ListIterate(I->Spec,tRec,next)) {
              if(onoff!=tRec->visible) {
                if(tRec->type==cExecObject) {
                  if(tRec->visible) {
                    tRec->in_scene = SceneObjectDel(G,tRec->obj);				
                    ExecutiveInvalidateSceneMembers(G);
                    tRec->visible=!tRec->visible;
                  } else {
                    if((!suppress_hidden)||(!hide_underscore)||(!tRec->is_hidden)) {
                      tRec->in_scene = SceneObjectAdd(G,tRec->obj);
                      ExecutiveInvalidateSceneMembers(G);
                      tRec->visible=!tRec->visible;
                    }
                  }
                } else if((tRec->type!=cExecSelection)||(!onoff)) /* hide all selections, but show all */
                  tRec->visible=!tRec->visible;
              }
            }
          }
          break;
        case cExecObject:
          /*
          if(rec->visible!=onoff) {
            if(rec->visible) {
              rec->in_scene = SceneObjectDel(G,rec->obj);				
              ExecutiveInvalidateSceneMembers(G);
            } else {
              rec->in_scene = SceneObjectAdd(G,rec->obj);
              ExecutiveInvalidateSceneMembers(G);
            }
            rec->visible=!rec->visible;
          }
          */
          if(onoff) { /* enable */
            ExecutiveSpecEnable(G,rec,parents,false);
          } else { /* disable */
            if(rec->visible) {
              if(rec->in_scene) 
                rec->in_scene = SceneObjectDel(G,rec->obj);				
              rec->visible = false;
              ExecutiveInvalidateSceneMembers(G);
            }
          }
          break;
        case cExecSelection:
          if(rec->visible!=onoff) {
            rec->visible=!rec->visible;
            if(rec->visible)
              if(SettingGetGlobal_b(G,cSetting_active_selections)) {
                ExecutiveHideSelections(G);
                rec->visible=true;
              }
            SceneInvalidate(G);
            SeqDirty(G);
          }
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
#else
  {
    SpecRec *tRec;

    if(strcmp(name,cKeywordAll)==0) {
      tRec=NULL;
      while(ListIterate(I->Spec,tRec,next)) {
        if(onoff!=tRec->visible) {
          if(tRec->type==cExecObject) {
            if(tRec->visible) {
              tRec->in_scene = SceneObjectDel(G,tRec->obj);				
              ExecutiveInvalidateSceneMembers(G);
            } else {
              tRec->in_scene = SceneObjectAdd(G,tRec->obj);
              ExecutiveInvalidateSceneMembers(G);
            }
          }
          if((tRec->type!=cExecSelection)||(!onoff)) /* hide all selections, but show all */
            tRec->visible=!tRec->visible;
        }
      }
    } else {
      tRec = ExecutiveFindSpec(G,name);
      if(tRec) {
        if(tRec->type==cExecObject) {
          if(tRec->visible!=onoff)
            {
              if(tRec->visible) {
                tRec->in_scene = SceneObjectDel(G,tRec->obj);				
                ExecutiveInvalidateSceneMembers(G);
              } else {
                tRec->in_scene = SceneObjectAdd(G,tRec->obj);
                ExecutiveInvalidateSceneMembers(G);
              }
              tRec->visible=!tRec->visible;
            }
        }
        else if(tRec->type==cExecSelection) {
          if(tRec->visible!=onoff) {
            tRec->visible=!tRec->visible;
            if(tRec->visible)
              if(SettingGetGlobal_b(G,cSetting_active_selections)) {
                ExecutiveHideSelections(G);
                tRec->visible=true;
              }
            SceneInvalidate(G);
            SeqDirty(G);
          }
        }
      }
    }
  }
#endif
  PRINTFD(G,FB_Executive)
    " ExecutiveSetObjVisib: leaving...\n"
    ENDFD;
  return 1;
}

/*========================================================================*/
void ExecutiveFullScreen(PyMOLGlobals *G,int flag)
{
  if(flag<0)
    flag = ! SettingGetGlobal_b(G,cSetting_full_screen);
#ifndef _PYMOL_NO_GLUT
 {
   register CExecutive *I = G->Executive;
   if(G->HaveGUI && G->ValidContext) {
     if(!SettingGet(G,cSetting_full_screen)) {
       I->oldPX = p_glutGet(P_GLUT_WINDOW_X) 
#ifdef FREEGLUT
         - p_glutGet(P_GLUT_WINDOW_BORDER_WIDTH)
#endif
         ;
       I->oldPY = p_glutGet(P_GLUT_WINDOW_Y) 
#ifdef FREEGLUT
         - p_glutGet(P_GLUT_WINDOW_HEADER_HEIGHT)
#endif
         ;
       I->oldWidth = p_glutGet(P_GLUT_WINDOW_WIDTH);
       I->oldHeight = p_glutGet(P_GLUT_WINDOW_HEIGHT);
       I->sizeFlag = true;
     }
     
     SettingSet(G,cSetting_full_screen,(float)flag);
     if(flag) {
       p_glutFullScreen();
     } else {
       if(I->sizeFlag) {
         p_glutPositionWindow(I->oldPX,I->oldPY);
         p_glutReshapeWindow(I->oldWidth,I->oldHeight);
       } else {
#ifndef _PYMOL_NO_MAIN
         MainRepositionWindowDefault(G);
#endif
       }
     }
   }
 }
#endif
  SettingSet(G,cSetting_full_screen,(float)flag);
  if(flag) {
     PyMOL_NeedReshape(G->PyMOL,1,0,0,0,0); /* zoom full-screen */
  } else {
     PyMOL_NeedReshape(G->PyMOL,0,0,0,0,0); /* return to non-zoomed size */
  }
  SceneChanged(G);
}
/*========================================================================*/
void ExecutiveSetAllVisib(PyMOLGlobals *G,int state)
{
  ObjectMoleculeOpRec op;
  ObjectMolecule *obj;
  int rep;
  int sele;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  PRINTFD(G,FB_Executive)
    " ExecutiveSetAllVisib: entered.\n"
    ENDFD;


  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject)
      {
        switch(rec->obj->type) {
        case cObjectMolecule:
          obj=(ObjectMolecule*)rec->obj;
          sele = SelectorIndexByName(G,obj->Obj.Name);
          for(rep=0;rep<cRepCnt;rep++) 
            rec->repOn[rep]=state;
          ObjectMoleculeOpRecInit(&op);

          op.code=OMOP_VISI;
          op.i1=-1;
          op.i2=state;
          ObjectMoleculeSeleOp(obj,sele,&op);
          op.code=OMOP_INVA;
          op.i1=-1;
          op.i2=cRepInvVisib;
          ObjectMoleculeSeleOp(obj,sele,&op);				
          break;
        default:
          for(rep=0;rep<cRepCnt;rep++) {
            ObjectSetRepVis(rec->obj,rep,state);
            if(rec->obj->fInvalidate)
              rec->obj->fInvalidate(rec->obj,rep,cRepInvVisib,state);
          }
          SceneInvalidate(G);
          break;
        }
      }
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveSetAllVisib: leaving...\n"
    ENDFD;

}
/*========================================================================*/
int ExecutiveToggleRepVisib(PyMOLGlobals *G,char *name,int rep)
{
  int ok =true;
  int sele;
  int handled = false;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;

  PRINTFD(G,FB_Executive)
    " ExecutiveToggleRepVisib: entered.\n"
    ENDFD;

  tRec = ExecutiveFindSpec(G,name);
  if((!tRec)&&(!strcmp(name,cKeywordAll))) {
    ExecutiveToggleAllRepVisib(G,rep);
  }
  if(tRec) {
    if(tRec->type==cExecObject) 
      switch(tRec->obj->type) {
      case cObjectMolecule: /* do nothing -- yet */
        break;
      default: /* otherwise, toggle the representation on/off */
        if(rep>=0) {
          ObjectToggleRepVis(tRec->obj,rep);
          if(tRec->obj->fInvalidate)
            tRec->obj->fInvalidate(tRec->obj,rep,cRepInvVisib,0);
        } 
        handled = true;
        SceneChanged(G);
        break;
      }
    if(!handled)
      switch(tRec->type) {
      case cExecSelection:
      case cExecObject:
        sele=SelectorIndexByName(G,name);
        if(sele>=0) {
          ObjectMoleculeOpRecInit(&op);

          op.code=OMOP_CheckVis;
          op.i1=rep;
          op.i2=false;
          ExecutiveObjMolSeleOp(G,sele,&op);
          op.i2 = !op.i2;

          if(tRec->type==cExecObject)
            ObjectSetRepVis(tRec->obj,rep,op.i2);

          op.code=OMOP_VISI;
          op.i1=rep;
          ExecutiveObjMolSeleOp(G,sele,&op);
          op.code=OMOP_INVA;
          op.i2=cRepInvVisib;
          ExecutiveObjMolSeleOp(G,sele,&op);
        }
        break;
      }
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveToggleRepVisib: leaving...\n"
    ENDFD;
  return (ok);
}

/*========================================================================*/
void ExecutiveSetRepVisib(PyMOLGlobals *G,char *name,int rep,int state)
{
  PRINTFD(G,FB_Executive)
    " ExecutiveSetRepVisib: entered.\n"
    ENDFD;

#if 1
  {
    register CExecutive *I = G->Executive;
    CTracker *I_Tracker= I->Tracker;
    SpecRec *rec = NULL;
    int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
    int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec) {
        /* per-atom */

        switch(rec->type) {
        case cExecSelection:
        case cExecObject: 
          {
            int sele=SelectorIndexByName(G,rec->name);
            if(sele>=0) {
              ObjectMoleculeOpRec op;
              ObjectMoleculeOpRecInit(&op);
              op.code=OMOP_VISI;
              op.i1=rep;
              op.i2=state;
              ExecutiveObjMolSeleOp(G,sele,&op);
              op.code=OMOP_INVA;
              op.i2=cRepInvVisib;
              ExecutiveObjMolSeleOp(G,sele,&op);
            }
          }
          break;
        }

        /* per-object/name */
        
        switch(rec->type) {
        case cExecSelection: 
          if(rec->name[0]!='_') {
            int a;
            /* remember visibility information for real selections */
            if(rep>=0) {
              rec->repOn[rep]=state;
            } else {
              for(a=0;a<cRepCnt;a++)
                rec->repOn[a]=state; 
            }
          }
          break;
        case cExecObject:
          if(rep>=0) {
            ObjectSetRepVis(rec->obj,rep,state);
            if(rec->obj->fInvalidate)
              rec->obj->fInvalidate(rec->obj,rep,cRepInvVisib,0);
          } else {
            int a;
            for(a=0;a<cRepCnt;a++) {
              rec->repOn[a]=state; 
              ObjectSetRepVis(rec->obj,a,state);
              if(rec->obj->fInvalidate)
                rec->obj->fInvalidate(rec->obj,a,cRepInvVisib,0);
            }
          }
          SceneChanged(G);
          break;
        case cExecAll:
          ExecutiveSetAllRepVisib(G,rep,state);
          break;
        }
      }
    }
    TrackerDelList(I_Tracker, list_id);
    TrackerDelIter(I_Tracker, iter_id);
  }
#else
 {
  int sele;
  int a;
  int handled = false;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;


  tRec = ExecutiveFindSpec(G,name);
  if((!tRec)&&(!strcmp(name,cKeywordAll))) {
    ExecutiveSetAllRepVisib(G,rep,state);
  }
  if(tRec) {
    if(name[0]!='_') {
      /* remember visibility information for real selections */
      if(rep>=0) {
        tRec->repOn[rep]=state;
      } else {
        for(a=0;a<cRepCnt;a++)
          tRec->repOn[a]=state; 
      }
    }
    if(tRec->type==cExecObject) 
      switch(tRec->obj->type) {
      default:
        if(rep>=0) {
          ObjectSetRepVis(tRec->obj,rep,state);
          if(tRec->obj->fInvalidate)
            tRec->obj->fInvalidate(tRec->obj,rep,cRepInvVisib,0);
        } else {
          for(a=0;a<cRepCnt;a++) {
            tRec->repOn[a]=state; 
            ObjectSetRepVis(tRec->obj,a,state);
            if(tRec->obj->fInvalidate)
              tRec->obj->fInvalidate(tRec->obj,a,cRepInvVisib,0);
          }
        }
        SceneChanged(G);
        break;
      }
    if(!handled)
      switch(tRec->type) {
      case cExecSelection:
      case cExecObject:
        sele=SelectorIndexByName(G,name);
        if(sele>=0) {
          ObjectMoleculeOpRecInit(&op);

          op.code=OMOP_VISI;
          op.i1=rep;
          op.i2=state;
          ExecutiveObjMolSeleOp(G,sele,&op);
          op.code=OMOP_INVA;
          op.i2=cRepInvVisib;
          ExecutiveObjMolSeleOp(G,sele,&op);
        }
        break;
      }
  }
 }
#endif

  PRINTFD(G,FB_Executive)
    " ExecutiveSetRepVisib: leaving...\n"
    ENDFD;

}

/*========================================================================*/
int ExecutiveSetOnOffBySele(PyMOLGlobals *G,char *name,int onoff)
{
  int sele;
  SpecRec *tRec;
  ObjectMoleculeOpRec op;

  tRec = ExecutiveFindSpec(G,name);
  if((!tRec)&&(!strcmp(name,cKeywordAll))) {
    ExecutiveSetObjVisib(G,name,onoff,false);
  }
  if(tRec) {
    sele=SelectorIndexByName(G,name);
    if(sele>=0) {
      ObjectMoleculeOpRecInit(&op);
      
      op.code=OMOP_OnOff;
      op.i1=onoff;
      ExecutiveObjMolSeleOp(G,sele,&op);
    }
  }
  return 1;
}


/*========================================================================*/
static void ExecutiveSetAllRepVisib(PyMOLGlobals *G,int rep,int state)
{
  ObjectMoleculeOpRec op;
  ObjectMolecule *obj;
  int sele;
  int a;

  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  PRINTFD(G,FB_Executive)
    " ExecutiveSetAllRepVisib: entered.\n"
    ENDFD;
  while(ListIterate(I->Spec,rec,next)) {
	 if(rec->type==cExecObject)
		{
        if(rec->name[0]!='_') {
          /* remember visibility information for real selections */
          if(rep>=0) {
            rec->repOn[rep]=state;
          } else {
            for(a=0;a<cRepCnt;a++)
              rec->repOn[a]=state; 
          }
        }   
        if(rec->type==cExecObject) {
          switch(rec->obj->type) {
          case cObjectMolecule:
            if(rep>=0) {
              rec->repOn[rep]=state;
            } else {
              for(a=0;a<cRepCnt;a++)
                rec->repOn[a]=state;
            }
            obj=(ObjectMolecule*)rec->obj;
            sele = SelectorIndexByName(G,obj->Obj.Name);
            ObjectMoleculeOpRecInit(&op);

            op.code=OMOP_VISI;
            op.i1=rep;
            op.i2=state;
            ObjectMoleculeSeleOp(obj,sele,&op);
            op.code=OMOP_INVA;
            op.i2=cRepInvVisib;
            ObjectMoleculeSeleOp(obj,sele,&op);				
            break;
          default:
            if(rep>=0) {
              ObjectSetRepVis(rec->obj,rep,state);
              if(rec->obj->fInvalidate)
                rec->obj->fInvalidate(rec->obj,rep,cRepInvVisib,state);
            } else {
              for(a=0;a<cRepCnt;a++) {
                ObjectSetRepVis(rec->obj,a,state);
                if(rec->obj->fInvalidate)
                  rec->obj->fInvalidate(rec->obj,rep,cRepInvVisib,state);
              }
            }
            SceneInvalidate(G);
            break;
          }
        }
		}
  }
  PRINTFD(G,FB_Executive)
    " ExecutiveSetAllRepVisib: leaving...\n"
    ENDFD;

}
/*========================================================================*/
static void ExecutiveToggleAllRepVisib(PyMOLGlobals *G,int rep)
{
  ObjectMoleculeOpRec op;
  int toggle_state;
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;

  op.code=OMOP_CheckVis;
  op.i1=rep;
  op.i2=false;
  ExecutiveObjMolSeleOp(G,cSelectionAll,&op);
  toggle_state = op.i2;
  while(ListIterate(I->Spec,rec,next)) {
	 if(rec->type==cExecObject) {
      switch(rec->obj->type) {
        case cObjectMolecule:
        break;
        default:
        if(rec->repOn[rep])
          toggle_state = true;
        break;
      }
    }
  }

  ExecutiveSetAllRepVisib(G,rep,!toggle_state);
}
/*========================================================================*/
void ExecutiveInvalidateRep(PyMOLGlobals *G,char *name,int rep,int level)
{
  register CExecutive *I = G->Executive;
  ObjectMoleculeOpRec op;
  SpecRec *rec = NULL;
  if((!name)||(!name[0])) 
    name = cKeywordAll;
#if 1
  {
      CTracker *I_Tracker= I->Tracker;
      int list_id = ExecutiveGetNamesListFromPattern(G,name,true,true);
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
        if(rec) {
          switch(rec->type) {
          case cExecSelection: 
          case cExecObject: 
            {
              int sele=SelectorIndexByName(G,rec->name);
              if(sele>=0) {
                ObjectMoleculeOpRecInit(&op);
                op.code = OMOP_INVA;
                op.i1=rep;
                op.i2=level;
                ExecutiveObjMolSeleOp(G,sele,&op);
              } else if(rec->obj->fInvalidate) {
                rec->obj->fInvalidate(rec->obj,rep,level,-1);
              }
            }
            break;
          case cExecAll:
            rec = NULL;
            while(ListIterate(I->Spec,rec,next)) {
              if(rec->type==cExecObject) {
                if(rec->obj->fInvalidate) {
                  rec->obj->fInvalidate(rec->obj,rep,level,-1);
                  SceneInvalidate(G);
                }
              }
            }
            break;
          }
        }
      }
      TrackerDelList(I_Tracker, list_id);
      TrackerDelIter(I_Tracker, iter_id);
  }
#else
  int sele = -1;
  int all_flag=false;
  PRINTFD(G,FB_Executive)
    "ExecInvRep-Debug: %s %d %d\n",name,rep,level
    ENDFD;
  if(WordMatch(G,cKeywordAll,name,true)<0) {
    all_flag=true;
  }
  if(all_flag) {
    while(ListIterate(I->Spec,rec,next))
      if(rec->type==cExecObject) {
        if(rec->obj->fInvalidate) {
          rec->obj->fInvalidate(rec->obj,rep,level,cRepAll);
          SceneInvalidate(G);
        }
      }
  }
  sele=SelectorIndexByName(G,name);
  if(sele>=0) {
    ObjectMoleculeOpRecInit(&op);
    op.code = OMOP_INVA;
    op.i1=rep;
    op.i2=level;
    ExecutiveObjMolSeleOp(G,sele,&op);
  }
#endif
}

int ExecutiveCheckGroupMembership(PyMOLGlobals *G,int list_id,CObject *obj)
{
  register CExecutive *I = G->Executive;
  int result = false;
  CTracker *I_Tracker= I->Tracker;
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  if(iter_id) {
    SpecRec *rec = NULL;
    while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
      if(rec && (rec->type==cExecObject) && (rec->obj == obj)) {
        result = true;
        break;
      }
    }
    TrackerDelIter(I_Tracker, iter_id);  
  }
  return result;
}

int ExecutiveGetExpandedGroupListFromPattern(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  int result = 0;
  CWordMatcher *matcher;
  CWordMatchOptions options;
  CTracker *I_Tracker= I->Tracker;
  char *wildcard = SettingGetGlobal_s(G,cSetting_wildcard);
  int iter_id = TrackerNewIter(I_Tracker, 0, I->all_names_list_id);
  int cand_id;
  SpecRec *rec;
      
  WordMatchOptionsConfigNameList(&options, 
                              *wildcard,
                              SettingGetGlobal_b(G,cSetting_ignore_case));
  matcher = WordMatcherNew(G, name, &options, false);
  if(matcher) {
    if(iter_id) {
      while( (cand_id = TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec)) ) {
        if(rec && !(rec->type==cExecAll)) {
          if(WordMatcherMatchAlpha(matcher,rec->name)) {
            if((rec->type==cExecObject)&&(rec->obj->type==cObjectGroup)) {
              if(!result)  
                result = TrackerNewList(I_Tracker, NULL);
              if(result) {
                TrackerLink(I_Tracker, cand_id, result, 1);
              }
            }
          }
        }
      }
    }
  } else if( (rec = ExecutiveFindSpec(G,name)) ) { /* only one name in list */
    if((rec->type==cExecObject) && (rec->obj->type==cObjectGroup)) {
      result = TrackerNewList(I_Tracker, NULL);
      TrackerLink(I_Tracker, rec->cand_id, result, 1);
    }
  } else if( (rec = ExecutiveUnambiguousNameMatch(G,name))) {
    if((rec->type==cExecObject) && (rec->obj->type==cObjectGroup)) {
      result = TrackerNewList(I_Tracker, NULL);
      TrackerLink(I_Tracker, rec->cand_id, result, 1);
    }
  }
  if(matcher) WordMatcherFree(matcher);
  if(iter_id) TrackerDelIter(I->Tracker, iter_id);
  if(result) 
    ExecutiveExpandGroupsInList(G,result,cExecExpandGroups);
  return result;
}

int ExecutiveGetExpandedGroupList(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  int result = 0;
  int list_id = 0;

  SpecRec *rec = ExecutiveFindSpec(G,name);
  ExecutiveUpdateGroups(G,false);
  if(rec && (rec->type == cExecObject) && (rec->obj->type == cObjectGroup )) {
    list_id = rec->group_member_list_id;
  }
  if(list_id) {
    result = TrackerNewListCopy(I->Tracker,list_id,NULL);
    ExecutiveExpandGroupsInList(G,result,cExecExpandGroups);
  }
  return result;
}

void ExecutiveFreeGroupList(PyMOLGlobals *G,int list_id)
{
  register CExecutive *I = G->Executive;
  TrackerDelList(I->Tracker, list_id);  
}


/*========================================================================*/
CObject *ExecutiveFindObjectByName(PyMOLGlobals *G,char *name)
{
  CObject *obj = NULL;
#if 1
  SpecRec *rec = ExecutiveFindSpec(G,name);
  if(rec && (rec->type == cExecObject)) {
    obj = rec->obj;
  }
#else
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  CObject *obj=NULL;
  while(ListIterate(I->Spec,rec,next)) { /* tisk tisk -- order N...
                                          * TODO use constant time lookup! */
    if(rec->type==cExecObject) {
      if(strcmp(rec->obj->Name,name)==0) {
        obj=rec->obj;
        break;
      }
    }
  }
#endif
  return(obj);
}
/*========================================================================*/
ObjectMap *ExecutiveFindObjectMapByName(PyMOLGlobals *G,char *name)
{
  CObject *obj;
  
  obj = ExecutiveFindObjectByName(G,name);
  if(obj)
    if(obj->type!=cObjectMap)
      obj=NULL;
  return((ObjectMap*)obj);
}

/*========================================================================*/
ObjectMolecule *ExecutiveFindObjectMoleculeByName(PyMOLGlobals *G,char *name)
{
  CObject *obj;
  
  obj = ExecutiveFindObjectByName(G,name);
  if(obj)
    if(obj->type!=cObjectMolecule)
      obj=NULL;
  return((ObjectMolecule*)obj);
}
/*========================================================================*/
Block *ExecutiveGetBlock(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  return(I->Block);
}
/*========================================================================*/
void ExecutiveSetControlsOff(PyMOLGlobals *G,char *name)
{
  SpecRec *rec;
  int a;
  rec = ExecutiveFindSpec(G,name);
  if(rec)
	 {
		for(a=0;a<cRepCnt;a++)
		  rec->repOn[a]=false;
	 }
}
/*========================================================================*/
void ExecutiveSymExp(PyMOLGlobals *G,char *name,
                     char *oname,char *s1,float cutoff,int segi,int quiet) /* TODO state */
{
  CObject *ob;
  ObjectMolecule *obj = NULL;
  ObjectMolecule *new_obj = NULL;
  ObjectMoleculeOpRec op;
  MapType *map;
  int x,y,z,a,b,c,i,j,h,k,l,n;
  CoordSet *cs,*os;
  int keepFlag,sele,tt[3];
  float *v1,*v2,m[16],tc[3],ts[3];
  OrthoLineType new_name;
  float auto_save;

  PRINTFD(G,FB_Executive)
    " ExecutiveSymExp: entered.\n"
    ENDFD;

  auto_save = SettingGet(G,cSetting_auto_zoom);
  SettingSet(G,cSetting_auto_zoom,0);
  sele=SelectorIndexByName(G,s1);
  ob = ExecutiveFindObjectByName(G,oname);
  if(ob->type==cObjectMolecule)
    obj=(ObjectMolecule*)ob;
  if(!(obj&&sele)) {
    ErrMessage(G,"ExecutiveSymExp","Invalid object");
  } else if(!obj->Symmetry) {
    ErrMessage(G,"ExecutiveSymExp","No symmetry loaded!");
  } else if(!obj->Symmetry->NSymMat) {
    ErrMessage(G,"ExecutiveSymExp","No symmetry matrices!");    
  } else {
    if(!quiet) {
      PRINTFB(G,FB_Executive,FB_Actions)
        " ExecutiveSymExp: Generating symmetry mates...\n"
        ENDFB(G);
    }
    ObjectMoleculeOpRecInit(&op);
	 op.code = OMOP_SUMC;
	 op.i1 =0;
    op.i2 =0;
    op.v1[0]= 0.0;
    op.v1[1]= 0.0;
    op.v1[2]= 0.0;
    ExecutiveObjMolSeleOp(G,sele,&op);
    tc[0]=op.v1[0];
    tc[1]=op.v1[1];
    tc[2]=op.v1[2];
    if(op.i1) {
      tc[0]/=op.i1;
      tc[1]/=op.i1;
      tc[2]/=op.i1;
    }
    transform33f3f(obj->Symmetry->Crystal->RealToFrac,tc,tc);

	 op.code = OMOP_VERT;
	 op.nvv1 =0;
    op.vv1 = VLAlloc(float,10000);
    ExecutiveObjMolSeleOp(G,sele,&op);
    
    if(!op.nvv1) {
      ErrMessage(G,"ExecutiveSymExp","No atoms indicated!");          
    } else {
      map=MapNew(G,-cutoff,op.vv1,op.nvv1,NULL);
      if(map) {
        MapSetupExpress(map);  
        /* go out no more than one lattice step in each direction */
        for(x=-1;x<2;x++)
          for(y=-1;y<2;y++)
            for(z=-1;z<2;z++)
              for(a=0;a<obj->Symmetry->NSymMat;a++) {
                new_obj = ObjectMoleculeCopy(obj);

                keepFlag=false;
                for(b=0;b<new_obj->NCSet;b++) 
                  if(new_obj->CSet[b]) {
                    cs = new_obj->CSet[b];
                    os = obj->CSet[b];
                    CoordSetRealToFrac(cs,obj->Symmetry->Crystal);
                    CoordSetTransform44f(cs,obj->Symmetry->SymMatVLA+(a*16));

                    CoordSetGetAverage(cs,ts);
                    identity44f(m);
                    /* compute the effective translation resulting
                       from application of the symmetry operator so
                       that we can shift it into the cell of the
                       target selection */
                    for(c=0;c<3;c++) { /* manual rounding - rint broken */
                      ts[c]=tc[c]-ts[c];
                      if(ts[c]<0)
                        ts[c]-=0.5;
                      else
                        ts[c]+=0.5;
                      tt[c]=(int)ts[c];
                    }
                    m[3] = (float)tt[0]+x;
                    m[7] = (float)tt[1]+y;
                    m[11] = (float)tt[2]+z;
                    CoordSetTransform44f(cs,m);
                    CoordSetFracToReal(cs,obj->Symmetry->Crystal);
                    if(!keepFlag) {
                      v2 = cs->Coord;
                      n=cs->NIndex;
                      while(n--) {
                        MapLocus(map,v2,&h,&k,&l);
                        i=*(MapEStart(map,h,k,l));
                        if(i) {
                          j=map->EList[i++];
                          while(j>=0) {
                            if(within3f(op.vv1+3*j,v2,cutoff)) {
                              keepFlag=true;
                              break;
                            }
                            j=map->EList[i++];
                          }
                        }
                        v2+=3;
                        if(keepFlag) break;
                      }
                    }
                    if(keepFlag) { /* make sure that we don't aren't simply duplicating the template coordinates */
                      keepFlag = false;
                      v1 = os->Coord;
                      v2 = cs->Coord;
                      n = cs->NIndex;
                      while(n--) {
                        if(diffsq3f(v1,v2)>R_SMALL8) {
                          keepFlag = true;
                          break;
                        }
                        v1++;
                        v2++;
                      }
                    }
                  }
                if(keepFlag) { /* we need to create new object */

                  /* TODO: should also transform the U tensor at this point...*/

                  sprintf(new_name,"%s%02d%02d%02d%02d",name,a,x,y,z);
                  ObjectSetName((CObject*)new_obj,new_name);
                  ExecutiveDelete(G,new_name);
                  ExecutiveManageObject(G,(CObject*)new_obj,-1,quiet);
                  SceneChanged(G);
                  if(segi==1) {
                    SegIdent seg;
                    if(a>35) {
                      seg[0] = 'a' + (a-36);
                    } else if(a>25) {
                      seg[0] = '0' + (a-26);
                    } else {
                      seg[0] = 'A' + a;
                    }
                    if(x>0) {
                      seg[1] = 'A' + x-1;
                    } else if(x<0) {
                      seg[1] = 'Z' + x+1;
                    } else {
                      seg[1] = '0';
                    }
                    if(y>0) {
                      seg[2] = 'A' + y-1;
                    } else if(y<0) {
                      seg[2] = 'Z' + y+1;
                    } else {
                      seg[2] = '0';
                    }
                    if(z>0) {
                      seg[3] = 'A' + z-1;
                    } else if(z<0) {
                      seg[3] = 'Z' + z+1;
                    } else {
                      seg[3] = '0';
                    }
                    seg[4] = 0;
                    {
                      int a;
                      AtomInfoType *ai = new_obj->AtomInfo;
                      for(a=0;a<new_obj->NAtom;a++) {
                        strcpy(ai->segi, seg);
                        ai++;
                      }

                    }
                  }
                } else {
                  ((CObject*)new_obj)->fFree((CObject*)new_obj);
                }
              }
        MapFree(map);
      }
    }
    VLAFreeP(op.vv1);
  }
  PRINTFD(G,FB_Executive)
     " ExecutiveSymExp: leaving...\n"
    ENDFD;
  SettingSet(G,cSetting_auto_zoom,auto_save);
}
static void ExecutivePurgeSpec(PyMOLGlobals *G,SpecRec *rec)
{
  register CExecutive *I = G->Executive;

  if(rec->group_name[0]) {
    /* cascade group members up to the surrounding group */
    SpecRec *rec2 = NULL;
    while(ListIterate(I->Spec,rec2,next)) {
      if((rec2->group == rec) ||
         WordMatch(G,rec->name,rec2->group_name,true)) {
        strcpy(rec2->group_name,rec->group_name);
      }
    }
  } else if((rec->type==cExecObject)&&(rec->obj->type == cObjectGroup)) {
    /* and/or delete their group membership */
    SpecRec *rec2 = NULL;
    while(ListIterate(I->Spec,rec2,next)) {
      if((rec2->group == rec) ||
         WordMatch(G,rec->name,rec2->group_name,true)) {
        rec2->group_name[0] = 0;
      }
    }
  }
  ExecutiveInvalidateGroups(G,false);
  ExecutiveInvalidatePanelList(G);
  switch(rec->type) {
  case cExecObject:
    if(I->LastEdited==rec->obj)
      I->LastEdited=NULL;
    if(rec->obj->type == cObjectMolecule)
      if(EditorIsAnActiveObject(G,(ObjectMolecule*)rec->obj))
        EditorInactivate(G);
    SeqChanged(G);
    if(rec->visible) { 
      SceneObjectDel(G,rec->obj);
      ExecutiveInvalidateSceneMembers(G);
    }
    ExecutiveDelKey(I,rec);
    SelectorDelete(G,rec->name);
    rec->obj->fFree(rec->obj);
    rec->obj=NULL;
    TrackerDelCand(I->Tracker,rec->cand_id);
    break;
  case cExecSelection:
    if(rec->visible) {
      SceneInvalidate(G);
      SeqDirty(G);
    }
    ExecutiveDelKey(I,rec);
    SelectorDelete(G,rec->name);
    TrackerDelCand(I->Tracker,rec->cand_id);
    break;
  }
}

/*========================================================================*/
void ExecutiveDelete(PyMOLGlobals *G,char *name)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
#if 1
  CTracker *I_Tracker= I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G,name,false,false);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
    if(rec) {
      switch(rec->type) {
      case cExecSelection: 
        ExecutivePurgeSpec(G,rec);
        ListDelete(I->Spec,rec,next,SpecRec); /* NOTE: order N in list length! TO FIX */
        break;
      case cExecObject: 
        if(rec->obj->type == cObjectGroup) {
          /* remove members of the group */
          ExecutiveGroup(G,rec->name,"",cExecutiveGroupPurge,true);
        }
        ExecutivePurgeSpec(G,rec);
        ListDelete(I->Spec,rec,next,SpecRec); /* NOTE: order N in list length! TO FIX */
        break;
      case cExecAll:
        rec = NULL;
        while(ListIterate(I->Spec, rec, next)) {
          switch(rec->type) {
          case cExecAll:
            break;
          default:
            ExecutivePurgeSpec(G,rec);
            ListDelete(I->Spec,rec,next,SpecRec);
            rec = NULL;
            break;
          }
        }
        SelectorDefragment(G);
        break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);
#else
  int all_flag = false;
  WordType name_copy; /* needed in case the passed string changes */

  if(WordMatch(G,name,cKeywordAll,true)<0) all_flag=true;
  strcpy(name_copy,name);
  while(ListIterate(I->Spec,rec,next)) {
    switch(rec->type) {
    case cExecObject:
      if(I->LastEdited==rec->obj)
        I->LastEdited=NULL;
      if(all_flag||(WordMatch(G,name_copy,rec->obj->Name,true)<0)) {
        
        if(rec->obj->type == cObjectMolecule)
          if(EditorIsAnActiveObject(G,(ObjectMolecule*)rec->obj))
            EditorInactivate(G);
        SeqChanged(G);
        if(rec->visible) {
          SceneObjectDel(G,rec->obj);
          ExecutiveInvalidateSceneMembers(G);
        }
        ExecutiveDelKey(I,rec);
        SelectorDelete(G,rec->name);
        rec->obj->fFree(rec->obj);
        rec->obj=NULL;
        TrackerDelCand(I->Tracker,rec->cand_id);
        ListDelete(I->Spec,rec,next,SpecRec);
        rec=NULL;
      }
      break;
    case cExecSelection:
      if(all_flag||(WordMatch(G,name_copy,rec->name,true)<0)) {
        
        if(all_flag||rec->visible)
          SceneInvalidate(G);
        SeqDirty(G);
        ExecutiveDelKey(I,rec);
        SelectorDelete(G,rec->name);
        TrackerDelCand(I->Tracker,rec->cand_id);
        ListDelete(I->Spec,rec,next,SpecRec);
        rec=NULL;
      }
      break;
    }
  }
  if(all_flag)
    SelectorDefragment(G);
#endif

}
/*========================================================================*/
void ExecutiveDump(PyMOLGlobals *G,char *fname,char *obj)
{
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;

  SceneUpdate(G,false);

  while(ListIterate(I->Spec,rec,next))
	 {
		if(rec->type==cExecObject)
		  {
			 if(strcmp(rec->obj->Name,obj)==0) 
				break;
		  }
	 }
  if(rec)
	 { 
      if(rec->obj->type==cObjectMesh) {
        ObjectMeshDump((ObjectMesh*)rec->obj,fname,0);
      } else if(rec->obj->type==cObjectSurface) {
        ObjectSurfaceDump((ObjectSurface*)rec->obj,fname,0);
      } else {
        ErrMessage(G,"ExecutiveDump","Invalid object type for this operation.");
      }
	 }
  else {
    ErrMessage(G,"ExecutiveDump","Object not found.");
  }
  
}
/*========================================================================*/
void ExecutiveMemoryDump(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  fprintf(stderr," Executive: %d candidate(s) %d list(s) %d link(s).\n",
          TrackerGetNCand(I->Tracker),
          TrackerGetNList(I->Tracker),
          TrackerGetNLink(I->Tracker));
}

void ExecutiveDoZoom(PyMOLGlobals *G,CObject *obj,int is_new, int zoom,int quiet)
{
  if(zoom) {/* -1 = use setting, 0 = never, 1 = zoom new, 
               2 = zoom always, 3 = zoom current, 4 = zoom all, 5= first object */
    if(zoom<0) {
      zoom = SettingGetGlobal_i(G,cSetting_auto_zoom);
      if(zoom<0) {
        zoom = 1;
      }
    } 
    switch(zoom) {
    case 1: /* zoom when new */
      if(is_new) 
        ExecutiveWindowZoom(G,obj->Name,0.0,-1,0,0,quiet); /* (all states) */
      break;
    case 2: /* zoom always */
      ExecutiveWindowZoom(G,obj->Name,0.0,-1,0,0,quiet); /* (all states) */
      break;
    case 3: /* always zoom current state */
      ExecutiveWindowZoom(G,obj->Name,0.0,ObjectGetCurrentState(obj,false),0,0,quiet); /* (all states) */
      break;
    case 4: /* zoom all objects */
      ExecutiveWindowZoom(G,cKeywordAll,0.0,-1,0,0,quiet);        
      break;
    case 5: /* zoom first object only */
      if(count_objects(G,true)==1)
        ExecutiveWindowZoom(G,obj->Name,0.0,-1,0,0,quiet); /* (all states) */ 
      break;
    }
  }
}
static void ExecutiveDoAutoGroup(PyMOLGlobals *G,SpecRec *rec)
{
  register CExecutive *I = G->Executive;
  int auto_mode = SettingGet(G,cSetting_group_auto_mode);
  if(auto_mode&&(rec->name[0]!='_')) {
    char *period = rec->name + strlen(rec->name);
    SpecRec *found_group = NULL;
    WordType seek_group_name;
    UtilNCopy(seek_group_name,rec->name,sizeof(WordType));
    
    while((period>rec->name)&&(!found_group)) {
      period--;
      if(*period=='.') {
        seek_group_name[period-rec->name]=0;
        {
          SpecRec *group_rec = NULL;
          while(ListIterate(I->Spec,group_rec,next)) {
            if((group_rec->type==cExecObject)&&(group_rec->obj->type==cObjectGroup)) {
              if(WordMatchExact(G,group_rec->name,seek_group_name,true)) {
                found_group = group_rec;
                strcpy(rec->group_name,seek_group_name);
                break;
              }
            }
          }
        }
        
        if((!found_group)&&(auto_mode==2)) {
          CObject *obj = (CObject*)ObjectGroupNew(G);
          if(obj) {
            ObjectSetName(obj,seek_group_name);
            strcpy(rec->group_name,seek_group_name);
            ExecutiveManageObject(G,obj,false,true);
            ExecutiveInvalidateGroups(G,false);
            break;
          }
        }
      }
    }
    if(found_group)
      ExecutiveInvalidateGroups(G,false);
  }
}

/*========================================================================*/
void ExecutiveManageObject(PyMOLGlobals *G,CObject *obj,int zoom,int quiet)
{
  int a;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;
  int exists=false;
  
  if(SettingGet(G,cSetting_auto_hide_selections))
    ExecutiveHideSelections(G);
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->obj==obj) {
      exists = true;
    }
  }
  if(!exists) {
    while(ListIterate(I->Spec,rec,next)) {
      if(rec->type==cExecObject) {
        if(strcmp(rec->obj->Name,obj->Name)==0) 
          break;
      }
    }
    if(rec) { /* another object of this type already exists */
      /* purge it */
      SceneObjectDel(G,rec->obj);
      ExecutiveInvalidateSceneMembers(G);
      rec->obj->fFree(rec->obj);
      rec->obj=NULL;
    } else {
      if(!quiet)
        if(obj->Name[0]!='_') { /* suppress internal objects */
          PRINTFB(G,FB_Executive,FB_Actions)
            " Executive: object \"%s\" created.\n",obj->Name 
            ENDFB(G);
        }
    }
    if(!rec)
      ListElemCalloc(G,rec,SpecRec);

    if(WordMatch(G,cKeywordAll,obj->Name,true)<0) {
      PRINTFB(G,FB_Executive,FB_Warnings) 
        " Executive: object name \"%s\" is illegal -- renamed to 'all_'.\n",obj->Name
        ENDFB(G);
      strcat(obj->Name,"_"); /* don't allow object named "all" */
    }
    if(SelectorNameIsKeyword(G, obj->Name)) {
      PRINTFB(G,FB_Executive,FB_Warnings) 
        " Executive-Warning: name \"%s\" collides with a selection language keyword.\n",obj->Name
        ENDFB(G);
    }
    strcpy(rec->name,obj->Name);
    rec->type=cExecObject;
    rec->next=NULL;
    rec->obj=obj;
    if(rec->obj->type==cObjectMap) {
      rec->visible=0;
    } else {
      rec->visible=1;
      /*      SceneObjectAdd(G,obj);*/
    }
    for(a=0;a<cRepCnt;a++)
      rec->repOn[a]=false;
    if(rec->obj->type==cObjectMolecule)
      rec->repOn[cRepLine]=true;

    rec->cand_id = TrackerNewCand(I->Tracker,(TrackerRef*)rec);

    TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id,1);
    TrackerLink(I->Tracker, rec->cand_id, I->all_obj_list_id,1);
    ListAppend(I->Spec,rec,next,SpecRec);
    ExecutiveAddKey(I,rec);
    ExecutiveInvalidatePanelList(G);

    if(rec->visible) {
      rec->in_scene = SceneObjectAdd(G,obj);
      ExecutiveInvalidateSceneMembers(G);
    }
    ExecutiveDoAutoGroup(G,rec);
  }

  if(obj->type==cObjectMolecule) {
	 ExecutiveUpdateObjectSelection(G,obj);
  }

  if(SettingGet(G,cSetting_auto_dss)) {
    if(obj->type==cObjectMolecule) {
      ObjectMolecule *objMol = (ObjectMolecule*)obj;
      if(objMol->NCSet==1) {
        ExecutiveAssignSS(G,obj->Name,0,NULL,1,1);
      }
    }
  }
  
  if(obj->fGetNFrame) {
    int n_state = obj->fGetNFrame(obj);
    int defer_limit = SettingGetGlobal_i(G,cSetting_auto_defer_builds);
    if( (defer_limit>=0) && (n_state >= defer_limit)) {
      int defer_builds = SettingGetGlobal_b(G,cSetting_defer_builds_mode);
      if(!defer_builds)
        SettingSetGlobal_b(G,cSetting_defer_builds_mode, 1);
    }
  }

  ExecutiveDoZoom(G,obj,!exists,zoom,true);

  SeqChanged(G);
}
/*========================================================================*/
void ExecutiveManageSelection(PyMOLGlobals *G,char *name)
{

  int a;
  SpecRec *rec = NULL;
  register CExecutive *I = G->Executive;
  int hide_all  = SettingGetGlobal_b(G,cSetting_active_selections);
  if(name[0]=='_')
    hide_all = false; /* hidden selections don't change active selection */
  while(ListIterate(I->Spec,rec,next))
    {
      if(rec->type==cExecSelection) {
        if(strcmp(rec->name,name)==0) 
          break;
        if(hide_all)
          rec->visible=false;
      }
    }
  if(rec&&hide_all)
    while(ListIterate(I->Spec,rec,next))
      if(rec->type==cExecSelection)
        rec->visible=false;

  if(!rec) {
    ListElemCalloc(G,rec,SpecRec);
    strcpy(rec->name,name);
    rec->type=cExecSelection;
    rec->next=NULL;
    rec->sele_color=-1;
    rec->visible=false;

    rec->cand_id = TrackerNewCand(I->Tracker,(TrackerRef*)rec);
    TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id,1);
    TrackerLink(I->Tracker, rec->cand_id, I->all_sel_list_id,1);
    ListAppend(I->Spec,rec,next,SpecRec);
    ExecutiveAddKey(I,rec);
    ExecutiveInvalidatePanelList(G);
  }
  if(rec) {
    for(a=0;a<cRepCnt;a++)
      rec->repOn[a]=false;
    if(name[0]!='_') {
      if(SettingGet(G,cSetting_auto_hide_selections))
        ExecutiveHideSelections(G);
      if(SettingGet(G,cSetting_auto_show_selections)) {
        rec->visible=true;
      }
    }
    if(rec->visible) SceneInvalidate(G);
    ExecutiveDoAutoGroup(G,rec);
  }
  SeqDirty(G);
}
/*========================================================================*/
static int ExecutiveClick(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CExecutive *I = G->Executive;
  int n,a;
  SpecRec *rec = NULL;
  PanelRec *panel = NULL;
  int t,xx;
  int pass = false;
  int skip;
  int ExecLineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);
  int hide_underscore = SettingGetGlobal_b(G,cSetting_hide_underscore_names);
  int op_cnt = get_op_cnt(G);

  if(y<I->HowFarDown) {
    if(SettingGetGlobal_b(G,cSetting_internal_gui_mode)==1) 
      return SceneDeferClick(SceneGetBlock(G),button,x,y,mod);
  }
  n=((I->Block->rect.top-y)-(ExecTopMargin+ExecClickMargin))/ExecLineHeight;
  a=n;
  xx = (x-I->Block->rect.left);
  if(I->ScrollBarActive) {
    if((x-I->Block->rect.left)<(ExecScrollBarWidth+ExecScrollBarMargin+ExecToggleMargin)) {
      pass = 1;
      ScrollBarDoClick(I->ScrollBar,button,x,y,mod);      
    }
    xx -= (ExecScrollBarWidth+ExecScrollBarMargin);
  } 
  skip = I->NSkip;
  if(!pass) {
    I->RecoverPressed = NULL;
    /* while(ListIterate(I->Spec,rec,next)) {*/
    while(ListIterate(I->Panel,panel,next)) {
      rec = panel->spec;

      if((rec->name[0]!='_')||(!hide_underscore)) {
        if(skip) {
          skip--;
        } else {
          if(!a) {
            t = ((I->Block->rect.right-ExecRightMargin)-x-1)/ExecToggleWidth;
            if(t<op_cnt) {
              int my = I->Block->rect.top-(ExecTopMargin + n*ExecLineHeight)-3;
              int mx = I->Block->rect.right-(ExecRightMargin + t*ExecToggleWidth);
              t = (op_cnt-t)-1;
              switch(t) {
              case 0: /* action */
                switch(rec->type) {
                case cExecAll:
                  MenuActivate(G,mx,my,x,y,false,"all_action",rec->name);
                  break;
                case cExecSelection:
                  MenuActivate(G,mx,my,x,y,false,"sele_action",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G,mx,my,x,y,false,"mol_action",rec->obj->Name);
                    break;
                  case cObjectMap:
                    MenuActivate(G,mx,my,x,y,false,"map_action",rec->obj->Name);
                    break;
                  case cObjectSurface:
                    MenuActivate(G,mx,my,x,y,false,"surface_action",rec->obj->Name);
                    break;
                  case cObjectMesh:
                    MenuActivate(G,mx,my,x,y,false,"mesh_action",rec->obj->Name);
                    break;
                  case cObjectMeasurement:
                  case cObjectCGO:
                  case cObjectCallback:
                  case cObjectAlignment:
                    MenuActivate(G,mx,my,x,y,false,"simple_action",rec->obj->Name);
                    break;
                  case cObjectSlice:
                    MenuActivate(G,mx,my,x,y,false,"slice_action",rec->obj->Name);
                    break;
                  case cObjectGadget:
                    MenuActivate(G,mx,my,x,y,false,"ramp_action",rec->obj->Name);
                    break;
                  }
                  break;
                }
                break;
              case 1:
                switch(rec->type) {
                case cExecAll:
                  MenuActivate(G,mx,my,x,y,false,"mol_show",cKeywordAll);
                  break;
                case cExecSelection:
                  MenuActivate(G,mx,my,x,y,false,"mol_show",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G,mx,my,x,y,false,"mol_show",rec->obj->Name);
                    break;
                  case cObjectCGO:
                  case cObjectAlignment:
                    MenuActivate(G,mx,my,x,y,false,"cgo_show",rec->obj->Name);
                    break;
                  case cObjectMeasurement:
                    MenuActivate(G,mx,my,x,y,false,"measurement_show",rec->obj->Name);
                    break;
                  case cObjectMap:
                    MenuActivate(G,mx,my,x,y,false,"map_show",rec->obj->Name);
                    break;
                  case cObjectMesh:
                    MenuActivate(G,mx,my,x,y,false,"mesh_show",rec->obj->Name);
                    break;
                  case cObjectSurface:
                    MenuActivate(G,mx,my,x,y,false,"surface_show",rec->obj->Name);
                    break;
                  case cObjectSlice:
                    MenuActivate(G,mx,my,x,y,false,"slice_show",rec->obj->Name);
                    break;
                    
                  }
                  break;
                }
                break;
              case 2:
                switch(rec->type) {
                case cExecAll:
                  MenuActivate(G,mx,my,x,y,false,"mol_hide",cKeywordAll);
                  break;
                case cExecSelection:
                  MenuActivate(G,mx,my,x,y,false,"mol_hide",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G,mx,my,x,y,false,"mol_hide",rec->obj->Name);
                    break;
                  case cObjectCGO:
                  case cObjectAlignment:
                    MenuActivate(G,mx,my,x,y,false,"cgo_hide",rec->obj->Name);
                    break;
                  case cObjectMeasurement:
                    MenuActivate(G,mx,my,x,y,false,"measurement_hide",rec->obj->Name);
                    break;
                  case cObjectMap:
                    MenuActivate(G,mx,my,x,y,false,"map_hide",rec->obj->Name);
                    break;
                  case cObjectMesh:
                    MenuActivate(G,mx,my,x,y,false,"mesh_hide",rec->obj->Name);
                    break;
                  case cObjectSurface:
                    MenuActivate(G,mx,my,x,y,false,"surface_hide",rec->obj->Name);
                    break;
                  case cObjectSlice:
                    MenuActivate(G,mx,my,x,y,false,"slice_hide",rec->obj->Name);
                    break;

                  }
                  break;
                }
                break;
              case 3:
                switch(rec->type) {
                case cExecAll:
                  MenuActivate(G,mx,my,x,y,false,"mol_labels","(all)");
                  break;
                case cExecSelection:
                  MenuActivate(G,mx,my,x,y,false,"mol_labels",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G,mx,my,x,y,false,"mol_labels",rec->obj->Name);
                    break;
                  case cObjectMeasurement:
                    break;
                  case cObjectMap:
                  case cObjectSurface:
                  case cObjectMesh:
                  case cObjectSlice:
                    break;
                  }
                  break;
                }
                break;
              case 4:
                switch(rec->type) {
                case cExecAll:
                case cExecSelection:
                  MenuActivate(G,mx,my,x,y,false,"mol_color",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G,mx,my,x,y,false,"mol_color",rec->obj->Name);
                    break;
                  case cObjectMeasurement:
                  case cObjectMap:
                  case cObjectSurface:
                  case cObjectCGO:
                  case cObjectMesh:
                    MenuActivate(G,mx,my,x,y,false,"general_color",rec->obj->Name);
                    break;
                  case cObjectSlice:
                    MenuActivate(G,mx,my,x,y,false,"slice_color",rec->obj->Name);
                    break;
                  }
                  break;
                }
                break;
              case 5:
                switch(rec->type) {
                case cExecAll:
                case cExecSelection:
                  MenuActivate(G,mx,my,x,y,false,"all_motion",rec->name);
                  break;
                case cExecObject:
                  switch(rec->obj->type) {
                  case cObjectGroup:
                  case cObjectMolecule:
                    MenuActivate(G,mx,my,x,y,false,"mol_motion",rec->obj->Name);
                    break;
                    /*
                      case cObjectMeasurement:
                      case cObjectMap:
                      case cObjectSurface:
                      case cObjectCGO:
                      case cObjectMesh:
                      MenuActivate(G,mx,my,x,y,false,"obj_motion",rec->obj->Name);
                      break;*/
                  }
                  break;
                }
                break;
              }
            } else { /* clicked in variable area */

              if(((panel->is_group)&&(((xx)-1)/8) > (panel->nest_level+1)) ||
                 ((!panel->is_group)&&(((xx)-1)/8) > panel->nest_level) ) {
                /* clicked on name */
                  
                rec->hilight=1;
                switch(button) {
                case P_GLUT_LEFT_BUTTON:
                  I->Pressed = n;
                  I->OldVisibility = rec->visible;
                  I->Over = n;
                  I->DragMode = 1;
                  I->ToggleMode=0;
                  I->LastChanged=NULL;
                  I->LastZoomed=NULL;
                  if(mod==(cOrthoSHIFT|cOrthoCTRL)) {
                    I->ToggleMode=2;
                    if(!rec->visible) {
                      ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod,false);
                    }
                    if(rec!=I->LastZoomed) 
                      ExecutiveWindowZoom(G,rec->name,0.0F,-1,false,-1.0F,true);
                    I->LastZoomed=rec;
                    I->LastChanged=rec;
                  } else if(mod&cOrthoSHIFT) {
                    ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod,false);
                    I->ToggleMode=1;
                  } else if(mod&cOrthoCTRL) {
                    I->ToggleMode=2;
                    if(!rec->visible) {
                      ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod,false);
                    }
                    I->LastChanged=rec;
                  }
                  I->PressedWhat = 1;
                  I->OverWhat = 1;
                  break;
                case P_GLUT_MIDDLE_BUTTON:
                  I->Pressed = n;
                  I->OldVisibility = rec->visible;
                  I->Over = n;
                  I->LastOver = I->Over;
                  I->DragMode = 3;
                  I->ToggleMode = 0;
                  I->LastChanged=rec;
                  I->LastZoomed=NULL;
                  if(mod&cOrthoCTRL) {
                    I->ToggleMode = 5;
                    ExecutiveWindowZoom(G,rec->name,0.0F,-1,false,-1.0F,true);
                    I->LastZoomed = rec;
                    if(mod&cOrthoSHIFT) { /* exclusive */
                      I->ToggleMode = 6;

                      ExecutiveSetObjVisib(G,cKeywordAll, false,false); /* need to log this */
                      if(!rec->visible)
                        ExecutiveSpecSetVisibility(G,rec,true,0,false);
                    }
                  } else {
                    I->ToggleMode = 4;
                    ExecutiveCenter(G,rec->name,-1,true,-1.0F,NULL,true);
                  }
                  if(!rec->visible) {
                    ExecutiveSpecSetVisibility(G,rec,!rec->visible,mod,false);
                    I->LastChanged=rec;
                  }
                  I->PressedWhat = 1;
                  I->OverWhat = 1;
                  break;
                case P_GLUT_RIGHT_BUTTON:
                  I->DragMode = 2; /* reorder */
                  I->Pressed = n;
                  I->Over = n;
                  I->LastOver = I->Over;
                  I->PressedWhat = 1;
                  I->OverWhat = 1;
                  break;
                }
                OrthoGrab(G,I->Block);
                OrthoDirty(G);
              } else if(panel->is_group) {
                /* clicked on group control */
                rec->hilight = 2;
                I->DragMode = 1;
                I->Pressed = n;
                I->Over = n;
                I->LastOver = I->Over;
                I->PressedWhat = 2;
                I->OverWhat = 2;

                OrthoGrab(G,I->Block);
                OrthoDirty(G);
              }
            }
          }
          a--;
        }
      }
    }
  }
  PyMOL_NeedRedisplay(G->PyMOL);
  return(1);
}
/*========================================================================*/
static void ExecutiveSpecEnable(PyMOLGlobals *G, SpecRec *rec, int parents, int log)
{
  if(log && SettingGetGlobal_b(G,cSetting_logging)) {
    OrthoLineType buffer = "";
    sprintf(buffer,"cmd.enable('%s',%d)",rec->obj->Name,parents);
    PLog(G,buffer,cPLog_pym);
  }
  
  if(!rec->visible) {
    rec->visible = true;
  }    
  if(!rec->in_scene) {
    rec->in_scene = SceneObjectAdd(G,rec->obj);
  }

  if(parents) {
    CExecutive *I = G->Executive;
    CTracker *I_Tracker= I->Tracker;
    int list_id = ExecutiveGetObjectParentList(G,rec);
    if(list_id) {
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec *parent_rec = NULL;
      
      while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&parent_rec) ) {
        if(rec) {
          switch(parent_rec->type) {
          case cExecObject:
            if(!parent_rec->in_scene) {
              parent_rec->in_scene = SceneObjectAdd(G,parent_rec->obj);                
            }
            if(!parent_rec->visible) {
              parent_rec->visible = true;
            }
          }
        }
      }
      TrackerDelIter(I_Tracker, iter_id);
    }
    TrackerDelList(I_Tracker, list_id);
  }
  ExecutiveInvalidateSceneMembers(G);
}

static void ExecutiveSpecSetVisibility(PyMOLGlobals *G,SpecRec *rec,
                                      int new_vis,int mod,int parents)
{
  OrthoLineType buffer = "";
  int logging = SettingGet(G,cSetting_logging);
  if(rec->type==cExecObject) {
    if(rec->visible&&!new_vis) {
      if(logging)
        sprintf(buffer,"cmd.disable('%s')",rec->obj->Name);
      SceneObjectDel(G,rec->obj);			
      ExecutiveInvalidateSceneMembers(G);
      rec->visible=new_vis;
    }
    else if((!rec->visible)&&new_vis) {
      ExecutiveSpecEnable(G,rec,parents,logging);
      /*
        if(logging) 
        sprintf(buffer,"cmd.enable('%s')",rec->obj->Name);
        rec->in_scene = SceneObjectAdd(G,rec->obj);
        rec->visible=new_vis;
        ExecutiveInvalidateSceneMembers(G);
      */
    }
    SceneChanged(G);
    if(logging&&buffer[0]) {
      PLog(G,buffer,cPLog_pym);
    }
  } else if(rec->type==cExecAll) {
    if(SettingGet(G,cSetting_logging)) {
      if(rec->visible)
        sprintf(buffer,"cmd.disable('all')");
      else
        sprintf(buffer,"cmd.enable('all')");
      PLog(G,buffer,cPLog_pym);
    }
    ExecutiveSetObjVisib(G,cKeywordAll,!rec->visible,false);
    } else if(rec->type==cExecSelection) {
      if(mod&cOrthoCTRL) {
        /*        SettingSet(G,cSetting_selection_overlay,
                  (float)(!((int)SettingGet(G,cSetting_selection_overlay))));
                  if(SettingGet(G,cSetting_logging)) {
                  sprintf(buffer,"cmd.set('selection_overlay',%d)",
                  (int)SettingGet(G,cSetting_selection_overlay));
                  PLog(G,buffer,cPLog_pym);
                  }
        */
        sprintf(buffer,"cmd.enable('%s')",rec->name);
        PLog(G,buffer,cPLog_pym);
        rec->visible=true; 
      } else if(mod&cOrthoSHIFT) {
        if(rec->sele_color<7)
          rec->sele_color=15;
        else {
          rec->sele_color--;
          if(rec->sele_color<7)
            rec->sele_color=15;
        }
        /* NO COMMAND EQUIVALENT FOR THIS FUNCTION YET */
        rec->visible=true;
      } else {
        
        if(rec->visible&&!new_vis) {
          if(SettingGet(G,cSetting_logging)) 
            sprintf(buffer,"cmd.disable('%s')",rec->name);
        }
        else if((!rec->visible)&&new_vis) {
          sprintf(buffer,"cmd.enable('%s')",rec->name);
        }
        if(new_vis && SettingGetGlobal_b(G,cSetting_active_selections)) {
          ExecutiveHideSelections(G);
        }
        if(SettingGet(G,cSetting_logging)) {
          PLog(G,buffer,cPLog_pym);
        }
        rec->visible=new_vis;
      }
      SceneChanged(G);
    }
}

static int ExecutiveRelease(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CExecutive *I = G->Executive;
  int n;  
  SpecRec *rec = NULL;
  PanelRec *panel = NULL;
  int pass = false;
  int skip;
  int xx;
  int ExecLineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);
  int hide_underscore = SettingGetGlobal_b(G,cSetting_hide_underscore_names);
  if(y<I->HowFarDown) {
    if(SettingGetGlobal_b(G,cSetting_internal_gui_mode)==1) 
      return SceneDeferRelease(SceneGetBlock(G),button,x,y,mod);
  }

  n=((I->Block->rect.top-y)-(ExecTopMargin+ExecClickMargin))/ExecLineHeight;

  xx = (x-I->Block->rect.left);
  if(I->ScrollBarActive) {
    if((x-I->Block->rect.left)<(ExecScrollBarWidth+ExecScrollBarMargin+ExecToggleMargin)) {
      pass = 1;
      ScrollBarDoRelease(I->ScrollBar,button,x,y,mod);
      OrthoUngrab(G);
    }
    xx -= (ExecScrollBarWidth+ExecScrollBarMargin);
  }

  skip=I->NSkip;

  if(!pass)
    {
      ExecutiveDrag(block,x,y,mod); /* incorporate final changes in cursor position */
      switch(I->DragMode) {
      case 1:
        
        /*while(ListIterate(I->Spec,rec,next)) {*/
        while(ListIterate(I->Panel,panel,next)) {
          rec = panel->spec;

          if((rec->name[0]!='_')||(!hide_underscore))
            {
              if(skip) {
                skip--;
              } else {
                if((I->PressedWhat==1) &&
                   (((panel->is_group)&&((xx-1)/8) > (panel->nest_level+1)) ||
                   ((!panel->is_group)&&((xx-1)/8) > panel->nest_level)) ) {
                  /* over name */
                  if(rec->hilight==1) {
                    if(rec->type==cExecSelection) {
                      ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,0,false);                    
                    } else {
                      ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod,true);
                    }
                  }
                } else if((I->PressedWhat==2) && (panel->is_group)) {
                  if(rec->hilight==2) {
                    ObjectGroup *obj = (ObjectGroup*)rec->obj;
                    OrthoLineType buf2;
                    sprintf(buf2,"cmd.group(\"%s\",action='%s')\n",rec->obj->Name,
                            obj->OpenOrClosed ? "close" : "open" );
                    PLog(G,buf2,cPLog_no_flush);
                    ExecutiveGroup(G,rec->obj->Name,"",5,1);

                  }
                  /* over group control */
                }
              }
            }
        }
        break;
      case 2:
        if(I->ReorderFlag) {
          I->ReorderFlag=false;
          PLog(G,I->ReorderLog,cPLog_no_flush);
        }
        break;
      }
    }
  
  {
    SpecRec *rec=NULL;
    while(ListIterate(I->Spec,rec,next)) {
      rec->hilight=0;
    }
  }

  I->Over = -1;
  I->Pressed = -1;
  I->DragMode = 0;
  I->PressedWhat = 0;
  OrthoUngrab(G);
  PyMOL_NeedRedisplay(G->PyMOL);
  return(1);
}
/*========================================================================*/
static int ExecutiveDrag(Block *block,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  register CExecutive *I = G->Executive;
  int xx,t;
  int ExecLineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);
  int hide_underscore = SettingGetGlobal_b(G,cSetting_hide_underscore_names);
  int op_cnt = get_op_cnt(G);


  if(y<I->HowFarDown) {
    if(SettingGetGlobal_b(G,cSetting_internal_gui_mode)==1) 
      return SceneDeferDrag(SceneGetBlock(G),x,y,mod);
  }

  if(I->DragMode) {
    xx = (x-I->Block->rect.left);
    t = ((I->Block->rect.right-ExecRightMargin)-x)/ExecToggleWidth;
    if(I->ScrollBarActive) {
      xx -= (ExecScrollBarWidth+ExecScrollBarMargin);
    }
  
    {
      int row_offset;
      if((xx>=0)&&(t>=op_cnt)) {
        row_offset = ((I->Block->rect.top-y)-
                      (ExecTopMargin+ExecClickMargin))/ExecLineHeight;
        I->Over = row_offset;
      } else {
        I->Over = -1;
        row_offset = -1;
        {
          SpecRec *rec=NULL;
          while(ListIterate(I->Spec,rec,next))      
            rec->hilight=0;
        }
      }
      
      if(I->RecoverPressed) {
        SpecRec *rec = NULL;
        PanelRec *panel = NULL;
        int skip=I->NSkip;
        int row = 0;
        while(ListIterate(I->Panel,panel,next)) {
          rec = panel->spec;
          
          if((rec->name[0]!='_')||(!hide_underscore)) {
            if(skip) {
              skip--;
            } else {
              if(rec==I->RecoverPressed) {
                I->Pressed = row;
                I->RecoverPressed = false;
              }
              row++;
            }
          }
        }
      }

      if(I->PressedWhat == 2) {
        I->OverWhat = 0;
      }

      if(I->Over>=0) {
        SpecRec *rec = NULL;
        PanelRec *panel = NULL;
        int skip=I->NSkip;
        int row=0;
        switch(I->DragMode) {
        case 1:
        
          /*          while(ListIterate(I->Spec,rec,next)) {*/
          while(ListIterate(I->Panel,panel,next)) {
            rec = panel->spec;

            if((rec->name[0]!='_')||(!hide_underscore)) {
              if(skip) {
                skip--;
              } else {
                
                rec->hilight=0;
                if((I->PressedWhat==1) && /* name button */
                   (((row>=I->Over)&&(row<=I->Pressed))||
                   ((row>=I->Pressed)&&(row<=I->Over)))) {
                  
                  if(((panel->is_group)&&((xx-1)/8) > (panel->nest_level+1)) ||
                      ((!panel->is_group)&&((xx-1)/8) > panel->nest_level)) {
                    /* dragged over name */
                    
                    I->OverWhat = 1;
                    
                    switch(I->ToggleMode) {
                    case 0:
                      if(row||(row==I->Pressed))
                        rec->hilight=1;
                      break;
                    case 1:
                      if(row)
                        ExecutiveSpecSetVisibility(G,rec,!I->OldVisibility,mod,false);
                      break;
                    case 2:
                      if((row==I->Over)&&row) {
                        if(I->LastChanged!=rec) {
                          ExecutiveSpecSetVisibility(G,I->LastChanged,false,mod,false);
                        }
                        if(!rec->visible) {
                          ExecutiveSpecSetVisibility(G,rec,true,mod,false);
                          I->LastChanged=rec;
                        }
                        if((mod==(cOrthoSHIFT|cOrthoCTRL))) {
                          if(rec!=I->LastZoomed) 
                            ExecutiveWindowZoom(G,rec->name,0.0F,-1,false,-1.0F,true);
                          I->LastZoomed=rec;
                        }
                      }
                      break;
                    }
                  }
                } else if((row==I->Pressed)&&(I->PressedWhat==2)) {
                  if(!((panel->is_group)&&((xx-1)/8) > (panel->nest_level+1))) {
                    
                    /* on group control */
                    
                    I->OverWhat = 2;
                    if(I->PressedWhat == I->OverWhat) {
                      rec->hilight = 2;
                    }
                  } else {
                    I->OverWhat = 0;
                  }
                }
                row++;
              }
            }
          }
          break;
        case 2: /* right buttonBROKEN */
          {
            if((I->Over!=I->Pressed) && (I->LastOver!=I->Over)) {
              SpecRec *new_rec = NULL;
              SpecRec *mov_rec = NULL;
              PanelRec *panel = NULL;
              
              while(ListIterate(I->Panel,panel,next)) {
                rec = panel->spec;
                if((rec->name[0]!='_')||(!hide_underscore)) {
                  {
                    if(skip) {
                      skip--;
                    } else {
                      if(row==I->Pressed) {
                       mov_rec = rec;
                      }
                      if(row==I->Over) {
                        new_rec = rec;
                      } 
                      row++;
                    }
                  }
                }
              }
              {
                int group_flag = false;
                int order_flag = false;
                char *first = NULL, *second = NULL;
                int is_child = false;
                
                if(mov_rec && (!new_rec) && (I->Over>I->Pressed) && mov_rec->group) {
                  first = mov_rec->group->name;
                  second = mov_rec->name;
                  order_flag=true;
                  strcpy(mov_rec->group_name, mov_rec->group->group_name);
                  group_flag=true;
                } else if(mov_rec && new_rec) {
                  if(mov_rec == new_rec->group) {
                    /* do nothing when a group is dragged over one of its members */
                  } else {

                    if(I->Pressed<I->Over) { /* downward */
                      first = new_rec->name;
                      second = mov_rec->name;
                      order_flag=true;
                    } else { /* upward */
                      first = mov_rec->name;
                      second = new_rec->name;
                      order_flag=true;
                    }
                  
                    if(mov_rec->group == new_rec->group) { /* reordering within a group level */
                      if((new_rec->type==cExecObject)&&(new_rec->obj->type==cObjectGroup)) {
                        ObjectGroup *group = (ObjectGroup*)new_rec->obj;
                        if(group->OpenOrClosed && !is_child) {
                          /* put inside an open group */
                          strcpy(mov_rec->group_name, new_rec->name);                        
                          order_flag = false;
                        group_flag = true;
                        }
                      }
                    } else if((mov_rec->group != new_rec) &&
                              (new_rec->group != mov_rec)) {
                      
                      if((new_rec->type==cExecObject)&&(new_rec->obj->type==cObjectGroup)) {
                        ObjectGroup *group = (ObjectGroup*)new_rec->obj;
                        if(group->OpenOrClosed && !is_child) {
                          /* put inside group */
                          strcpy(mov_rec->group_name, new_rec->name);                        
                        } else if(new_rec->group_name) {
                          strcpy(mov_rec->group_name, new_rec->group_name);
                        } else {
                          mov_rec->group_name[0] = 0;
                        }
                        
                        /* WLD TEST */
                        
                        if(I->Pressed<I->Over) { /* downward */
                          first = new_rec->name;
                          second = mov_rec->name;
                          order_flag=true;
                        } else { /* upward */
                          first = mov_rec->name;
                          second = new_rec->name;
                          order_flag=true;
                        }

                      /* WLD END */

                    } else if(new_rec->group_name) 
                      strcpy(mov_rec->group_name, new_rec->group_name);
                      else
                        mov_rec->group_name[0] = 0;
                      group_flag = true;
                    }
                  }
                }
                
                if(group_flag) {
                  OrthoLineType buf;
                  if(mov_rec->group_name[0]) {
                    sprintf(buf,"group %s, %s\n",mov_rec->group_name,mov_rec->name);
                  } else {
                    sprintf(buf,"ungroup %s\n",mov_rec->name);
                  }
                  PLog(G,buf,cPLog_no_flush);                      
                  ExecutiveInvalidateGroups(G,false);
                  I->RecoverPressed = mov_rec;
                  I->Pressed = 0;
                }
                if(order_flag && first && second) {
                  OrthoLineType order_input;
                  sprintf(order_input,"%s %s",first,second);
                  ExecutiveOrder(G,order_input,false,0);
                  sprintf(I->ReorderLog,"cmd.order(\"%s\")\n",
                          order_input);
                  PLog(G,I->ReorderLog,cPLog_no_flush);
                I->RecoverPressed = mov_rec;
                I->Pressed = 0;
                }
              }
            }
          }
          break;
        case 3: /* middle button */
          /* while(ListIterate(I->Spec,rec,next)) {*/
          while(ListIterate(I->Panel,panel,next)) {
            rec = panel->spec;
            if((rec->name[0]!='_')||(!hide_underscore))
              {
                if(skip) {
                  skip--;
                } else {
                  rec->hilight=0;
                  if( ((row>=I->Over)&&(row<=I->Pressed))||
                      ((row>=I->Pressed)&&(row<=I->Over))) {
                    switch(I->ToggleMode) {
                    case 4: /* center and activate, while deactivating previous */
                      if((row==I->Over)&&row) {
                        if(I->LastChanged!=rec) {
                          ExecutiveSpecSetVisibility(G,I->LastChanged,false,mod,false);
                          ExecutiveCenter(G,rec->name,-1,true,-1.0F,NULL,true);
                          if(!rec->visible) 
                            ExecutiveSpecSetVisibility(G,rec,true,mod,false);
                          I->LastChanged = rec;
                        }
                        rec->hilight=0;
                      }
                      break;
                    case 5: /* zoom and activate, while deactivating previous */
                      if((row==I->Over)&&row) {
                        if(I->LastChanged!=rec) {
                          ExecutiveSpecSetVisibility(G,I->LastChanged,false,mod,false);
                          ExecutiveWindowZoom(G,rec->name,0.0F,-1,false,-1.0F,true);
                          if(!rec->visible) 
                            ExecutiveSpecSetVisibility(G,rec,true,mod,false);
                          I->LastChanged = rec;
                        }
                        rec->hilight=1;
                      }
                      break;
                    case 6: /* zoom and make only object enabled */
                      if((row==I->Over)&&row) {
                        if(rec!=I->LastZoomed) { 
                          ExecutiveSpecSetVisibility(G,I->LastZoomed,false,mod,false);
                          ExecutiveWindowZoom(G,rec->name,0.0F,-1,false,-1.0F,true);
                          I->LastZoomed=rec;
                          ExecutiveSpecSetVisibility(G,rec,true,0,false);
                        }
                        rec->hilight=1;
                      }
                    }
                  }
                  row++;
                }
              }
          }
          break;
        }
        I->LastOver = I->Over;
        
      } else if(I->LastChanged)
        ExecutiveSpecSetVisibility(G,I->LastChanged,false,mod,false);      
      OrthoDirty(G);
    }
  }
  return(1);
}

static void draw_button(int x2,int y2, int w, int h, float *light, float *dark, float *inside)
{
  glColor3fv(light);
  glBegin(GL_POLYGON);
  glVertex2i(x2,y2);
  glVertex2i(x2,y2+h);
  glVertex2i(x2+w,y2+h);
  glVertex2i(x2+w,y2);
  glEnd();
  
  glColor3fv(dark);
  glBegin(GL_POLYGON);
  glVertex2i(x2+1,y2);
  glVertex2i(x2+1,y2+h-1);
  glVertex2i(x2+w,y2+h-1);
  glVertex2i(x2+w,y2);
  glEnd();
  
  if(inside) {
    glColor3fv(inside);
    glBegin(GL_POLYGON);
    glVertex2i(x2+1,y2+1);
    glVertex2i(x2+1,y2+h-1);
    glVertex2i(x2+w-1,y2+h-1);
    glVertex2i(x2+w-1,y2+1);
    glEnd();
  } else { /* rainbow */
    glBegin(GL_POLYGON);
    glColor3f(1.0F,0.1F,0.1F);
    glVertex2i(x2+1,y2+1);
    glColor3f(0.1F,1.0F,0.1F);
    glVertex2i(x2+1,y2+h-1);
    glColor3f(1.0F,1.0F,0.1F);
    glVertex2i(x2+w-1,y2+h-1);
    glColor3f(0.1F,0.1F,1.0F);
    glVertex2i(x2+w-1,y2+1);
    glEnd();
  }

}

#ifndef _PYMOL_NOPY
static void draw_button_char(PyMOLGlobals *G,int x2,int y2,char ch)
{
  TextSetColor3f(G,0.0F,0.0F,0.0F);
  TextSetPos2i(G,x2+ExecToggleTextShift,y2);
  TextDrawChar(G,ch);
}
#endif


/*========================================================================*/
static void ExecutiveDraw(Block *block)
{
  PyMOLGlobals *G=block->G;
  int x,y,xx,x2,y2;
  char *c=NULL;
  float enabledColor[3] = { 0.5F, 0.5F, 0.5F };
  float cloakedColor[3] = { 0.35F, 0.35F, 0.35F };
  float pressedColor[3] = { 0.7F, 0.7F, 0.7F };
  float disabledColor[3] = { 0.25F, 0.25F, 0.25F };
  float lightEdge[3] = {0.6F, 0.6F, 0.6F };
  float darkEdge[3] = {0.35F, 0.35F, 0.35F };
  float captionColor[3] = {0.3F, 0.9F, 0.3F };
#ifndef _PYMOL_NOPY
  float toggleColor3[3] = { 0.6F, 0.6F, 0.8F };
#endif

  SpecRec *rec = NULL;
  PanelRec *panel = NULL;
  register CExecutive *I = G->Executive;
  int n_ent;
  int n_disp;
  int skip=0;
  int row = -1;
  int ExecLineHeight = SettingGetGlobal_i(G,cSetting_internal_gui_control_size);
  int text_lift = (ExecLineHeight/2)-5;
  int hide_underscore = SettingGetGlobal_b(G,cSetting_hide_underscore_names);
  int op_cnt = get_op_cnt(G);
  int full_names = SettingGetGlobal_b(G,cSetting_group_full_member_names);
  int arrows = SettingGetGlobal_b(G,cSetting_group_arrow_prefix);

  ExecutiveUpdatePanelList(G);
  if(G->HaveGUI && G->ValidContext && 
     ((block->rect.right-block->rect.left)>6) 
     && I->ValidPanel) {
    int max_char;
    int nChar;
    /* do we have enough structures to warrant a scroll bar? */
    
    n_ent = 0;
    while(ListIterate(I->Panel,panel,next)) {
      rec=panel->spec;
      if(rec && ((rec->name[0]!='_')||(!hide_underscore)))
        n_ent++;
    }

    n_disp = ((I->Block->rect.top-I->Block->rect.bottom)-(ExecTopMargin))/ExecLineHeight;
    if(n_disp<1) n_disp=1;
      
    if(n_ent>n_disp) {
      int bar_maxed = ScrollBarIsMaxed(I->ScrollBar);
      if(!I->ScrollBarActive) {
        ScrollBarSetLimits(I->ScrollBar,n_ent,n_disp);
        if(bar_maxed) {
          ScrollBarMaxOut(I->ScrollBar);
          I->NSkip = (int)ScrollBarGetValue(I->ScrollBar);
        } else {
          ScrollBarSetValue(I->ScrollBar,0);
          I->NSkip =0;
        }
      } else {
        ScrollBarSetLimits(I->ScrollBar,n_ent,n_disp);
        if(bar_maxed)
          ScrollBarMaxOut(I->ScrollBar);
        I->NSkip = (int)ScrollBarGetValue(I->ScrollBar);
      }
      I->ScrollBarActive = 1;

    } else {
      I->ScrollBarActive = 0;
      I->NSkip =0;
    }

    max_char = (((I->Block->rect.right-I->Block->rect.left)-(ExecLeftMargin+ExecRightMargin+4)) -
                     (op_cnt*ExecToggleWidth));
    if(I->ScrollBarActive) {
      max_char -= (ExecScrollBarMargin+ExecScrollBarWidth);
    }      
    max_char/=8;

    if(SettingGetGlobal_b(G,cSetting_internal_gui_mode)==0) {
      glColor3fv(I->Block->BackColor);
      BlockFill(I->Block);
    }

    if(I->ScrollBarActive) {
      ScrollBarSetBox(I->ScrollBar,I->Block->rect.top-ExecScrollBarMargin,
                      I->Block->rect.left+ExecScrollBarMargin,
                      I->Block->rect.bottom+2,
                      I->Block->rect.left+ExecScrollBarMargin+ExecScrollBarWidth);
      ScrollBarDoDraw(I->ScrollBar);
    }
    
    x = I->Block->rect.left+ExecLeftMargin;
    y = (I->Block->rect.top-ExecLineHeight)-ExecTopMargin;
    /*    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(cRepCnt+op_cnt);*/
#ifndef _PYMOL_NOPY
    xx = I->Block->rect.right-ExecRightMargin-ExecToggleWidth*(op_cnt);
#else
    xx = I->Block->rect.right-ExecRightMargin;
#endif

    if(I->ScrollBarActive) {
      x+=ExecScrollBarWidth+ExecScrollBarMargin;
    }
    skip=I->NSkip;
    /*    while(ListIterate(I->Spec,rec,next)) {*/
    while(ListIterate(I->Panel,panel,next)) {
      rec = panel->spec;
      if((rec->name[0]!='_')||(!hide_underscore))
        {
          if(skip) {
            skip--;
          } else {
            row++;
            x2=xx;
            y2=y;
            nChar = max_char;

            if((x-ExecToggleMargin)-(xx-ExecToggleMargin)>-10) {
              x2 = x+10;
            }
#ifndef _PYMOL_NOPY
            {
              int a;
              float toggleColor[3] = { 0.5F, 0.5F, 1.0F };
              float toggleColor2[3] = { 0.4F, 0.4F, 0.6F };
              float toggleDarkEdge[3] = { 0.3F, 0.3F, 0.5F};
              float toggleLightEdge[3] = { 0.7F, 0.7F, 0.9F};
              
              glColor3fv(toggleColor);
              for(a=0;a<op_cnt;a++)
                {
                  switch(a) {
                  case 0:
                    /*
                      glColor3fv(toggleColor);
                      glBegin(GL_POLYGON);
                      glVertex2i(x2,y2+(ExecToggleSize)/2);
                      glVertex2i(x2+(ExecToggleSize)/2,y2);
                      glVertex2i(x2+ExecToggleSize,y2+(ExecToggleSize)/2);
                      glVertex2i(x2+(ExecToggleSize)/2,y2+ExecToggleSize);
                      glEnd();
                    */
                    
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor);
                    
                    draw_button_char(G,x2,y2+text_lift,'A');
                    break;
                  case 1:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor3);
                    
                    draw_button_char(G,x2,y2+text_lift,'S');
                    break;
                  case 2:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor2);
                    draw_button_char(G,x2,y2+text_lift,'H');
                    break;
                  case 3:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor);
                    draw_button_char(G,x2,y2+text_lift,'L');
                    break;
                  case 4:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                NULL);
                    draw_button_char(G,x2,y2+text_lift,'C');
                    break;
                  case 5:
                    draw_button(x2,y2,ExecToggleSize,(ExecLineHeight-1),
                                toggleLightEdge,
                                toggleDarkEdge,
                                toggleColor3);
                    draw_button_char(G,x2,y2+text_lift,'M');
                    break;
                  }
                  x2+=ExecToggleWidth;
                }
            }
#endif

            {
              int x3 = x;
              int hidden_prefix = false;
              
              TextSetColor(G,I->Block->TextColor);
              TextSetPos2i(G,x3+2,y2+text_lift);
              
              if((rec->type==cExecObject)||
                 (rec->type==cExecAll)||
                 (rec->type==cExecSelection)) {

                y2=y;
                x2 = xx;
                if((x-ExecToggleMargin)-(xx-ExecToggleMargin)>-10) {
                  x2 = x+10;
                }
                x3+=panel->nest_level*8;
                TextSetPos2i(G,x3+2,y2+text_lift);
                nChar-=panel->nest_level;
                {
                  int but_width = (x2-x3)-1;

                  
                  if(panel->is_group) {

                    if((rec->hilight==2)&&(I->Over==I->Pressed)) {
                      draw_button(x3,y2,15,(ExecLineHeight-1),lightEdge,darkEdge,pressedColor);
                    } else if(panel->is_open) {
                      draw_button(x3,y2,15,(ExecLineHeight-1),lightEdge,darkEdge,disabledColor);
                    } else {
                      draw_button(x3,y2,15,(ExecLineHeight-1),lightEdge,darkEdge,disabledColor);
                    }
                    
#define cControlBoxSize 17
#define cControlLeftMargin 8
#define cControlTopMargin 2
#define cControlSpacing 2
#define cControlInnerMargin 4
#define cControlSpread 6
#define cControlSize 160

                    TextSetPos2i(G,x3+4,y2+text_lift);
                    if(panel->is_open) {
                      TextDrawChar(G,'-');
                    } else {
                      TextDrawChar(G,'+');
                    }

                    but_width-=16;
                    x3+=16;
                    nChar-=2;

                    TextSetPos2i(G,x3+2,y2+text_lift);
                  }
                  
                  if((rec->hilight==1)||((row==I->Over)&&(I->OverWhat==1))) {
                    draw_button(x3,y2,but_width,(ExecLineHeight-1),lightEdge,darkEdge,pressedColor);
                  } else if(rec->visible) {
                    int enabled = true;
                    SpecRec *group_rec = rec->group;
                    while(enabled && group_rec ) { 
                      if(!group_rec->visible)
                        enabled=false;
                      else
                        group_rec = group_rec->group;
                    }

                    if(enabled) {
                      draw_button(x3,y2,but_width,(ExecLineHeight-1),lightEdge,darkEdge,enabledColor);
                    } else {
                      draw_button(x3,y2,but_width,(ExecLineHeight-1),lightEdge,darkEdge,cloakedColor);
                    }
                  } else {
                    draw_button(x3,y2,but_width,(ExecLineHeight-1),lightEdge,darkEdge,disabledColor);
                  }
                }
                
                TextSetColor(G,I->Block->TextColor);
                
                c=rec->name;

                if(!full_names) {
                  if(rec->group) { /* if prefix matches group name, then truncate */
                    char *p=c,* q=rec->group->name;
                    while((*p==*q)&&(*q)) {
                      p++;
                      q++;
                    }
                    if((*p)&&(!*q)&&(*p=='.')) {
                      hidden_prefix = true;
                      c=p;
                    }
                  }
                }

                if(rec->type==cExecSelection)
                  if((nChar--)>0) {
                    TextDrawChar(G,'(');
                  }
              }
              
              if(c) {
                if(hidden_prefix) {
                  if(arrows&&((nChar--)>0)) {
                    TextDrawChar(G,'^');
                    TextSetPos2i(G,x3+2,y2+text_lift);
                    TextDrawChar(G,'|');
                  }
                }
                  
                while(*c) {
                  if((nChar--)>0)
                    TextDrawChar(G,*(c++));
                  else
                    break;
                }
              }
              
              if(rec->type==cExecSelection)
                {
                  if((nChar--)>0) {
                    TextDrawChar(G,')');
                  }
                  
                  c=rec->name;
                }
              
              if(rec->type==cExecObject) {
                if(rec->obj->fGetCaption)
                  c = rec->obj->fGetCaption(rec->obj);
                if(c && c[0] && nChar>1 && strcmp(c,rec->obj->Name)!=0) {
                  TextSetColor(G,captionColor);
                  TextSetPos2i(G,x+2+8*(max_char-nChar),y2+text_lift);
                  if((nChar--)>0)
                    TextDrawChar(G,' ');
                  while(*c) 
                    if((nChar--)>0) 
                      TextDrawChar(G,*(c++));
                    else
                      break;
                }
              }
            }
            
            y-=ExecLineHeight;
            if(y<(I->Block->rect.bottom))
              break;
          }
        }
  }
    I->HowFarDown = y;
  }
}
/*========================================================================*/
int ExecutiveIterateObject(PyMOLGlobals *G,CObject **obj,void **hidden)
{
  int result;
  register CExecutive *I = G->Executive;
  int flag=false;
  SpecRec **rec=(SpecRec**)hidden;
  while(!flag) {
    result = (ListIterate(I->Spec,(*rec),next))!=NULL;
    if(!(*rec))
      flag=true;
    else if((*rec)->type==cExecObject)
      flag=true;
  }
  if(*rec)
	 (*obj)=(*rec)->obj;
  else
	 (*obj)=NULL;
  return(result);
}
/*========================================================================*/
int ExecutiveIterateObjectMolecule(PyMOLGlobals *G,ObjectMolecule **obj,void **hidden)
{
  int result;
  register CExecutive *I = G->Executive;
  int flag=false;
  SpecRec **rec=(SpecRec**)hidden;
  while(!flag)
	 {
		result = (ListIterate(I->Spec,(*rec),next))!=NULL;
		if(!(*rec))
		  flag=true;
		else if((*rec)->type==cExecObject)
        if((*rec)->obj->type==cObjectMolecule)
          flag=true;
	 }
  if(*rec)
	 (*obj)=(ObjectMolecule*)(*rec)->obj;
  else
	 (*obj)=NULL;
  return(result);
}
/*========================================================================*/
static void ExecutiveReshape(Block *block,int width,int height)
{
  PyMOLGlobals *G=block->G;
  register CExecutive *I = G->Executive;

  BlockReshape(block,width,height);

  I->Width = block->rect.right-block->rect.left+1;
  I->Height = block->rect.top-block->rect.bottom+1;
  
}

/*========================================================================*/
int ExecutiveReinitialize(PyMOLGlobals *G,int what,char *pattern)
{ 
  register CExecutive *I = G->Executive;
  int ok=true;
#ifndef _PYMOL_NOPY      
  int blocked = false;
#endif
  /* reinitialize PyMOL */
  if(what==2)
    pattern = NULL;

  if(pattern&&(!pattern[0])) pattern=NULL;
  if(!pattern) {
    
    switch(what) {
    case 0: /* everything */
      ExecutiveDelete(G,cKeywordAll);
      ColorReset(G);
      SettingInitGlobal(G,false,false,true);
      MovieReset(G);
      EditorInactivate(G);
      ControlRock(G,0);

#ifndef _PYMOL_NOPY      
      blocked = PAutoBlock(G);
      PRunStringInstance(G,"cmd.view('*','clear')");
      PRunStringInstance(G,"cmd.scene('*','clear')");
      WizardSet(G,NULL,false);
      PAutoUnblock(G,blocked);
#endif

      SculptCachePurge(G);
      SceneReinitialize(G);
      SelectorReinit(G);
      SeqChanged(G);
      break;
    case 1: /* settings */
      SettingInitGlobal(G,false,false,true);
      ExecutiveRebuildAll(G);
      break;
    case 2: /* store_defaults */
      SettingStoreDefault(G);
      break;
    case 3: /* original_settings */
      SettingInitGlobal(G,false,false,false);
      ExecutiveRebuildAll(G);
      break;
    case 4: /* purge_defaults */
      SettingPurgeDefault(G);
      break;
    }
  } else {
    {
      CTracker *I_Tracker= I->Tracker;
      int list_id = ExecutiveGetNamesListFromPattern(G,pattern,true,true);
      int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
      SpecRec *rec;
      
      while( TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef**)&rec) ) {
        if(rec) {
          switch(rec->type) {
          case cExecObject:
            switch(what) {
            case 0:
            case 1:
              if(rec->obj->Setting) {
                ObjectPurgeSettings(rec->obj);
                if(rec->obj->fInvalidate)
                  rec->obj->fInvalidate(rec->obj,cRepAll,cRepInvAll,-1);
                SceneInvalidate(G);
                SeqChanged(G);
              }
              break;
            }
          }
        }
      }
      TrackerDelList(I_Tracker, list_id);
      TrackerDelIter(I_Tracker, iter_id);
    }
    
    /* to do */
  }
  return(ok);
}
/*========================================================================*/
int ExecutiveInit(PyMOLGlobals *G)
{
  register CExecutive *I=NULL;
  if( (I=(G->Executive=Calloc(CExecutive,1)))) {
    
    SpecRec *rec = NULL;
    int a;

    ListInit(I->Spec);
    I->Tracker = TrackerNew(G);
    I->all_names_list_id = TrackerNewList(I->Tracker, NULL);
    I->all_obj_list_id = TrackerNewList(I->Tracker,NULL);
    I->all_sel_list_id = TrackerNewList(I->Tracker,NULL);
    I->Block = OrthoNewBlock(G,NULL);  
    I->Block->fRelease = ExecutiveRelease;
    I->Block->fClick   = ExecutiveClick;
    I->Block->fDrag    = ExecutiveDrag;
    I->Block->fDraw    = ExecutiveDraw;
    I->Block->fReshape = ExecutiveReshape;
    I->Block->active = true;
    I->ScrollBarActive = 0;
    I->ScrollBar=ScrollBarNew(G,false);
    OrthoAttach(G,I->Block,cOrthoTool);
    I->RecoverPressed = NULL;
    I->Pressed = -1;
    I->Over = -1;
    I->LastEdited=NULL;
    I->ReorderFlag=false;
    I->NSkip=0;
    I->HowFarDown=0;
    I->DragMode = 0;
    I->sizeFlag=false;
    I->LastZoomed = NULL;
    I->LastChanged = NULL;
    I->ValidGroups = false;
    I->ValidSceneMembers = false;

    ListInit(I->Panel);
    I->ValidPanel = false;

    I->Lex = OVLexicon_New(G->Context->heap);
    I->Key = OVOneToOne_New(G->Context->heap);

    /* create "all" entry */

    ListElemCalloc(G,rec,SpecRec);
  
    strcpy(rec->name,cKeywordAll);
    rec->type=cExecAll;
    rec->visible=true;
    rec->next=NULL;
    for(a=0;a<cRepCnt;a++)
      rec->repOn[a]=false;
    rec->cand_id = TrackerNewCand(I->Tracker,(TrackerRef*)rec);
    TrackerLink(I->Tracker, rec->cand_id, I->all_names_list_id,1);
    ListAppend(I->Spec,rec,next,SpecRec);
    ExecutiveAddKey(I,rec);    

    return 1;
  }
  else return 0;

}
/*========================================================================*/
void ExecutiveFree(PyMOLGlobals *G)
{
  register CExecutive *I = G->Executive;
  SpecRec *rec=NULL;
  while(ListIterate(I->Spec,rec,next)) {
    if(rec->type==cExecObject)
      rec->obj->fFree(rec->obj);
  }
  ListFree(I->Spec,next,SpecRec);
  ListFree(I->Panel,next,PanelRec);
  if(I->Tracker)
    TrackerFree(I->Tracker);
  if(I->ScrollBar)
    ScrollBarFree(I->ScrollBar);
  OrthoFreeBlock(G,I->Block);
  I->Block=NULL;
  OVLexicon_DEL_AUTO_NULL(I->Lex);
  OVOneToOne_DEL_AUTO_NULL(I->Key);

  FreeP(G->Executive);
}


#ifdef _undefined

matrix checking code...

		double mt[3][3],mt2[3][3],pr[3][3],im[3][3],em[3][3];
	 printf("normalized matrix \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%12.3f ",evect[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 printf("tensor \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%12.3f ",mi[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	  for(a=0;a<3;a++) {
		 for(b=0;b<3;b++) {
			mt[a][b]=evect[a][b];
		 }
	  }

	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  mt2[a][b]=evect[b][a];
		}
	 }

	 matrix_multiply33d33d(mt,mt2,pr);
	 printf("self product 1 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 matrix_multiply33d33d(mt,mi,im);
	 matrix_multiply33d33d(im,mt2,pr);
	 printf("diagonal product 1 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");

	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  em[a][b]=0.0;
		}
		em[a][a]=egval[a];
	 }

	 matrix_multiply33d33d(mt2,em,im);
	 matrix_multiply33d33d(im,mt,pr);
	 printf("diagonal product 4 \n");
	 for(a=0;a<3;a++) {
		for(b=0;b<3;b++) {
		  printf("%8.3f ",pr[a][b]);
		}
		printf("\n");
	 }
	 printf("\n");
#endif


