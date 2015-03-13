
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
#include"os_python.h"
#include "os_std.h"
#include "MemoryDebug.h"
#include "Err.h"
#include "PlugIOManager.h"
#include "Selector.h"
#include "CoordSet.h"
#include "Feedback.h"
#include "Scene.h"
#include "Executive.h"
#include "AtomInfo.h"

#ifndef _PYMOL_VMD_PLUGINS
int PlugIOManagerInit(PyMOLGlobals * G)
{
  return 1;
}

int PlugIOManagerFree(PyMOLGlobals * G)
{
  return 1;
}

int PlugIOManagerRegister(PyMOLGlobals * G, void *ptr);
int PlugIOManagerRegister(PyMOLGlobals * G, void *ptr)
{
  return 1;
}

int PlugIOManagerLoadTraj(PyMOLGlobals * G, ObjectMolecule * obj,
                          const char *fname, int frame,
                          int interval, int average, int start,
                          int stop, int max, char *sele, int image,
                          float *shift, int quiet, const char *plugin_type)
{

  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " ObjectMolecule-Error: sorry, VMD Molfile Plugins not compiled into this build.\n"
    ENDFB(G);
  return 0;
}

ObjectMap *PlugIOManagerLoadVol(PyMOLGlobals * G, ObjectMap * obj,
                                const char *fname, int state, int quiet,
                                const char *plugin_type)
{
  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " ObjectMap-Error: sorry, VMD Molfile Plugins not compiled into this build.\n"
    ENDFB(G);
  return NULL;
}

ObjectMolecule *PlugIOManagerLoadMol(PyMOLGlobals * G, ObjectMolecule *origObj,
    const char *fname, int state, int quiet, const char *plugin_type)
{
  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " ObjectMolecule-Error: sorry, VMD Molfile Plugins not compiled into this build.\n"
    ENDFB(G);
  return 0;
}

CObject * PlugIOManagerLoad(PyMOLGlobals * G, CObject ** obj_ptr,
    const char *fname, int state, int quiet, const char *plugin_type)
{
  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " ObjectMolecule-Error: sorry, VMD Molfile Plugins not compiled into this build.\n"
    ENDFB(G);
  return 0;
}

#else

#include "molfile_plugin.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _CPlugIOManager {
  int NPlugin;
  molfile_plugin_t **PluginVLA;
};

int PlugIOManagerInitAll(PyMOLGlobals * G);     /* defined externally */

int PlugIOManagerInit(PyMOLGlobals * G)
{
  CPlugIOManager *I = NULL;
  if((I = (G->PlugIOManager = Calloc(CPlugIOManager, 1)))) {
    I->NPlugin = 0;
    I->PluginVLA = VLAlloc(molfile_plugin_t *, 10);
    return PlugIOManagerInitAll(G);
  } else
    return 0;
}

int PlugIOManagerFreeAll(void); /* defined externally */

int PlugIOManagerFree(PyMOLGlobals * G)
{
  CPlugIOManager *I = G->PlugIOManager;
  PlugIOManagerFreeAll();
  VLAFreeP(I->PluginVLA);
  FreeP(G->PlugIOManager);
  return 1;
}

int PlugIOManagerRegister(PyMOLGlobals * G, vmdplugin_t * header);

int PlugIOManagerRegister(PyMOLGlobals * G, vmdplugin_t * header)
{
  if(G && G->PlugIOManager) {
    if(!strcmp(header->type, MOLFILE_PLUGIN_TYPE)) {
      CPlugIOManager *I = G->PlugIOManager;
      VLACheck(I->PluginVLA, molfile_plugin_t *, I->NPlugin);
      I->PluginVLA[I->NPlugin] = (molfile_plugin_t *) header;
      I->NPlugin++;
      /*           printf("register %p %s\n",header,header->name); */
    }
    return VMDPLUGIN_SUCCESS;
  } else
    return VMDPLUGIN_ERROR;
}

static molfile_plugin_t * find_plugin(CPlugIOManager * I, const char * plugin_type) {
  for (int a = 0; a < I->NPlugin; a++)
    if(!strcmp(plugin_type, I->PluginVLA[a]->name))
      return I->PluginVLA[a];
  return NULL;
}

int PlugIOManagerLoadTraj(PyMOLGlobals * G, ObjectMolecule * obj,
                          const char *fname, int frame,
                          int interval, int average, int start,
                          int stop, int max, char *sele, int image,
                          float *shift, int quiet, const char *plugin_type)
{
  CPlugIOManager *I = G->PlugIOManager;
  molfile_plugin_t *plugin = NULL;

  ok_assert(1, I && obj);
  plugin = find_plugin(I, plugin_type);

  if(!plugin) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " PlugIOManager: unable to locate plugin '%s'\n", plugin_type ENDFB(G);
    return false;
  }

  if(plugin->read_next_timestep == NULL) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " PlugIOManager: not a trajectory plugin '%s'\n", plugin_type ENDFB(G);
    return false;
  }

  {
      int natoms;
      molfile_timestep_t timestep;
      void *file_handle;
      int zoom_flag = false;
      int cnt = 0;
      int icnt = interval;
      int n_avg = 0;
      int ncnt = 0;
      CoordSet *cs = obj->NCSet > 0 ? obj->CSet[0] : obj->CSTmpl ? obj->CSTmpl : NULL;

      timestep.coords = NULL;
      timestep.velocities = NULL;

      file_handle = plugin->open_file_read(fname, plugin_type, &natoms);

      if(!file_handle) {
        PRINTFB(G, FB_ObjectMolecule, FB_Errors)
          " ObjectMolecule: plugin '%s' cannot open '%s'.\n", plugin_type, fname ENDFB(G);
        return false;
      }

      if(natoms == -1) {
        natoms = obj->NAtom;
      } else if(natoms != obj->NAtom) {
	PRINTFB(G, FB_ObjectMolecule, FB_Errors)
          " ObjectMolecule: plugin '%s' cannot open file because the number "
          "of atoms in the object (%d) did not equal the number of atoms in "
          "the '%s' (%d) file.\n", plugin_type, obj->NAtom, plugin_type, natoms ENDFB(G);
        return false;
      }

      if(cs) {
        ok_assert(1, cs = CoordSetCopy(cs));
      } else {
        ok_assert(1, cs = CoordSetNew(G));
        ok_assert(1, cs->Coord = VLAlloc(float, 3 * natoms));

        cs->Obj = obj;
        cs->NIndex = natoms;
        cs->enumIndices();
      }

      timestep.coords = (float *) cs->Coord;

      {
	  /* read_next_timestep fills in &timestep for each iteration; we need
	   * to copy that out to a new CoordSet, each time. */
          while(!plugin->read_next_timestep(file_handle, natoms, &timestep)) {
            cnt++;
	    /* start at the 'start'-th frame; skip 'start' frames,
	     * and skip every interval/icnt frames */
            if(cnt >= start) {
              icnt--;
              if(icnt > 0) {
                PRINTFB(G, FB_ObjectMolecule, FB_Details)
                  " ObjectMolecule: skipping set %d...\n", cnt ENDFB(G);
              } else {
                icnt = interval;
                n_avg++;
              }
              if(icnt == interval) {
                if(n_avg < average) {
                  PRINTFB(G, FB_ObjectMolecule, FB_Details)
                    " ObjectMolecule: averaging set %d...\n", cnt ENDFB(G);
                } else {
                  /* compute average */
                  if(n_avg > 1) {
                    float *fp;
                    int i;
                    fp = cs->Coord;
                    for(i = 0; i < cs->NIndex; i++) {
                      *(fp++) /= n_avg;
                      *(fp++) /= n_avg;
                      *(fp++) /= n_avg;
                    }
                  }
                  /* add new coord set */
                  cs->invalidateRep(cRepAll, cRepInvRep);
                  if(frame < 0) frame = obj->NCSet;
                  if(!obj->NCSet) zoom_flag = true;

		  /* make sure we have room for 'frame' CoordSet*'s in obj->CSet */
		  /* TODO: TEST this function */
                  VLACheck(obj->CSet, CoordSet*, frame); /* was CoordSet* */
		  /* bump the object's state count */
                  if(obj->NCSet <= frame) obj->NCSet = frame + 1;
		  /* if there's data in this state's coordset, emtpy it */
                  if(obj->CSet[frame])
                    obj->CSet[frame]->fFree();
		  /* set this state's coordset to cs */
                  obj->CSet[frame] = cs;
                  ncnt++;
                  if(average < 2) {
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                      " ObjectMolecule: read set %d into state %d...\n", cnt, frame + 1
                      ENDFB(G);
                  } else {
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                      " ObjectMolecule: averaging set %d...\n", cnt ENDFB(G);
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                      " ObjectMolecule: average loaded into state %d...\n", frame + 1
                      ENDFB(G);
                  }

                  if((stop > 0 && cnt >= stop) || (max > 0 && ncnt >= max)) {
                    cs = NULL;
                    break;
                  }

                  frame++;
		  /* make a new cs */
                  cs = CoordSetCopy(cs);        /* otherwise, we need a place to put the next set */
                  timestep.coords = (float *) cs->Coord;
                  n_avg = 0;
                }
              }
            } else {
              PRINTFB(G, FB_ObjectMolecule, FB_Details)
                " ObjectMolecule: skipping set %d...\n", cnt ENDFB(G);
            }
          } /* end while */
        }
        plugin->close_file_read(file_handle);
        if(cs)
          cs->fFree();
        SceneChanged(G);
        SceneCountFrames(G);
        if(zoom_flag)
          if(SettingGetGlobal_i(G, cSetting_auto_zoom)) {
            ExecutiveWindowZoom(G, obj->Obj.Name, 0.0, -1, 0, 0, quiet);        /* auto zoom (all states) */
          }
  }
  return true;
ok_except1:
  return false;
}

ObjectMap *PlugIOManagerLoadVol(PyMOLGlobals * G, ObjectMap * obj,
                                const char *fname, int state, int quiet,
                                const char *plugin_type)
{
  CPlugIOManager *I = G->PlugIOManager;
  molfile_plugin_t *plugin = NULL;
  molfile_volumetric_t *metadata;
  int setsinfile = 0;
  int natoms;
  void *file_handle = NULL;
  float *datablock = NULL;

  ok_assert(1, I);
  plugin = find_plugin(I, plugin_type);

  if(!plugin) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " PlugIOManager: unable to locate plugin '%s'\n", plugin_type ENDFB(G);
    ok_raise(1);
  }

  if(plugin->read_volumetric_data == NULL || plugin->read_volumetric_metadata == NULL) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " PlugIOManager: not a map plugin '%s'\n", plugin_type ENDFB(G);
    ok_raise(1);
  }

  file_handle = plugin->open_file_read(fname, plugin_type, &natoms);

  if(!file_handle) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " PlugIOManager: plugin '%s' cannot open '%s'.\n", plugin_type, fname ENDFB(G);
    ok_raise(1);
  }

  if(plugin->read_volumetric_metadata(file_handle, &setsinfile, &metadata) != MOLFILE_SUCCESS) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " PlugIOManager: read_volumetric_metadata failed\n" ENDFB(G);
    ok_raise(1);
  }

  /* for now, just read in all sets of data in the file */

  for(int i = 0; i < setsinfile; i++) {

          const molfile_volumetric_t *v = metadata + i;
          size_t size = v->xsize * v->ysize * v->zsize;

          if(!size) {
            PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
              " PlugIOManagerLoadVol-Waring: 0 values, skipping set %d\n", i ENDFB(G);
            continue;
          }

          ok_assert(1, datablock = Alloc(float, size));

          if(plugin->read_volumetric_data(file_handle, i, datablock, NULL) != MOLFILE_SUCCESS) {
            PRINTFB(G, FB_ObjectMolecule, FB_Errors)
              " PlugIOManager: read_volumetric_data failed\n" ENDFB(G);
            ok_raise(1);
          }

          {
            ObjectMapState *ms = NULL;

            if(!obj)
              ok_assert(1, obj = ObjectMapNew(G));

            if(state < 0)
              state = obj->NState;
            if(obj->NState <= state) {
              VLACheck(obj->State, ObjectMapState, state);
              obj->NState = state + 1;
            }
            ms = &obj->State[state];
            ObjectMapStateInit(obj->Obj.G, ms);

            ms->FDim[0] = v->xsize;
            ms->FDim[1] = v->ysize;
            ms->FDim[2] = v->zsize;
            ms->FDim[3] = 3;

            ms->Grid = Alloc(float, 3);
            ms->Dim = Alloc(int, 3);
            ms->Origin = Calloc(float, 3);
            ms->Range = Alloc(float, 3);

            // set corners to a unit cube, and manage world space with the state matrix
            {
              double m44d[16];

              if(!ms->State.Matrix)
                ms->State.Matrix = Alloc(double, 16);

              // state matrix transformation
              identity44d(m44d);
              copy3(v->xaxis, m44d + 0);
              copy3(v->yaxis, m44d + 4);
              copy3(v->zaxis, m44d + 8);
              copy3(v->origin, m44d + 12);
              transpose44d44d(m44d, ms->State.Matrix);
            }

            // axis and corner stuff in a unit cube
            {
              // prime min+max
              zero3f(ms->ExtentMin);
              ones3f(ms->ExtentMax);
              ones3f(ms->Range);

              for(int a = 0; a < 3; a++) {
                int dimL1 = ms->FDim[a] - 1;

                // min+max indices
                ms->Min[a] = 0;
                ms->Max[a] = dimL1;

                ms->Grid[a] = 1.f / dimL1;

                // (redundant)
                ms->Dim[a] = ms->FDim[a];

                // corner enumeration
                for(int b = 0; b < 8; b++)
                  ms->Corner[3 * b + a] = (b >> a) & 0x1;
              }
            }

            // field
            ms->Field = IsosurfFieldAlloc(G, ms->FDim);
            ms->MapSource = cMapSourceVMDPlugin;
            ms->Field->save_points = false;     /* save points in RAM only, not session file */
            ms->Active = true;

            ObjectMapStateRegeneratePoints(ms);

            // copy data
            {
              int a, b, c;
              float *data_ptr = datablock;

              /* VMD plugins appear to use fast-x med-y slow-z ordering: "&datablock[z*xysize + y*xsize + x]" */

              for(c = 0; c < ms->FDim[2]; c++) {
                for(b = 0; b < ms->FDim[1]; b++) {
                  for(a = 0; a < ms->FDim[0]; a++) {
                    F3(ms->Field->data, a, b, c) = *(data_ptr++);
                  }
                }
              }

              PRINTFB(G, FB_ObjectMap, FB_Details)
                " ObjectMap: read %zu values\n", size ENDFB(G);
            }

          }
          FreeP(datablock);
  }
  if(obj) {
    ObjectMapUpdateExtents(obj);
    SceneChanged(G);
    SceneCountFrames(G);
  }

ok_except1:
  // close
  if (plugin && file_handle)
    plugin->close_file_read(file_handle);

  return obj;
}

static void atomicnumber2elem(char * dst, int protons) {
  const char * p = NULL;
  switch(protons) {
    case cAN_LP: p = "LP"; break;
    case cAN_H:  p = "H";  break;
    case cAN_He: p = "He"; break;
    case cAN_Li: p = "Li"; break;
    case cAN_Be: p = "Be"; break;
    case cAN_B:  p = "B";  break;
    case cAN_C:  p = "C";  break;
    case cAN_N:  p = "N";  break;
    case cAN_O:  p = "O";  break;
    case cAN_F:  p = "F";  break;
    case cAN_Ne: p = "Ne"; break;
    case cAN_Na: p = "Na"; break;
    case cAN_Mg: p = "Mg"; break;
    case cAN_Al: p = "Al"; break;
    case cAN_Si: p = "Si"; break;
    case cAN_P:  p = "P";  break;
    case cAN_S:  p = "S";  break;
    case cAN_Cl: p = "Cl"; break;
    case cAN_Ar: p = "Ar"; break;
    case cAN_K:  p = "K";  break;
    case cAN_Ca: p = "Ca"; break;
    case cAN_Ti: p = "Ti"; break;
    case cAN_Cr: p = "Cr"; break;
    case cAN_Mn: p = "Mn"; break;
    case cAN_Fe: p = "Fe"; break;
    case cAN_Co: p = "Co"; break;
    case cAN_Ni: p = "Ni"; break;
    case cAN_Cu: p = "Cu"; break;
    case cAN_Zn: p = "Zn"; break;
    case cAN_Ga: p = "Ga"; break;
    case cAN_Ge: p = "Ge"; break;
    case cAN_As: p = "As"; break;
    case cAN_Se: p = "Se"; break;
    case cAN_Br: p = "Br"; break;
    case cAN_Kr: p = "Kr"; break;
    case cAN_Rb: p = "Rb"; break;
    case cAN_Sr: p = "Sr"; break;
    case cAN_Pd: p = "Pd"; break;
    case cAN_Ag: p = "Ag"; break;
    case cAN_Cd: p = "Cd"; break;
    case cAN_In: p = "In"; break;
    case cAN_Sn: p = "Sn"; break;
    case cAN_Sb: p = "Sb"; break;
    case cAN_Te: p = "Te"; break;
    case cAN_I:  p = "I";  break;
    case cAN_Xe: p = "Xe"; break;
    case cAN_Cs: p = "Cs"; break;
    case cAN_Ba: p = "Ba"; break;
    case cAN_Ce: p = "Ce"; break;
    case cAN_Pt: p = "Pt"; break;
    case cAN_Au: p = "Au"; break;
    case cAN_Hg: p = "Hg"; break;
    case cAN_Tl: p = "Tl"; break;
    case cAN_Pb: p = "Pb"; break;
    case cAN_U:  p = "U";  break;
    default:     p = "X";
  }
  strncpy(dst, p, cElemNameLen);
}

static CSymmetry * SymmetryNewFromTimestep(PyMOLGlobals * G, molfile_timestep_t * ts)
{
  CSymmetry * symm = NULL;
  ok_assert(1,
      ts->A > 0.f && ts->B > 0.f && ts->C > 0.f &&
      ts->alpha > 0.f && ts->beta > 0.f && ts->gamma > 0.f);
  ok_assert(1, symm = SymmetryNew(G));
  symm->Crystal->Dim[0] = ts->A;
  symm->Crystal->Dim[1] = ts->B;
  symm->Crystal->Dim[2] = ts->C;
  symm->Crystal->Angle[0] = ts->alpha;
  symm->Crystal->Angle[1] = ts->beta;
  symm->Crystal->Angle[2] = ts->gamma;
  strcpy(symm->SpaceGroup, "P1");
  SymmetryAttemptGeneration(symm, false);
ok_except1:
  return symm;
}

ObjectMolecule *PlugIOManagerLoadMol(PyMOLGlobals * G, ObjectMolecule *origObj,
    const char *fname, int state, int quiet, const char *plugin_type)
{
  CPlugIOManager *manager = G->PlugIOManager;
  int natoms, nbonds = 0, *from, *to;
  int optflags = 0;
  float *order;
  void *file_handle = NULL;
  molfile_plugin_t * plugin = NULL;
  molfile_timestep_t timestep;
  molfile_atom_t * atoms = NULL;
  ObjectMolecule *I = NULL;
  CoordSet * cs = NULL;
  int *bondtype, nbondtypes;
  char **bondtypename;
  int auto_show = RepGetAutoShowMask(G);

  ok_assert(1, manager);
  plugin = find_plugin(manager, plugin_type);

  if (!plugin) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjectMolecule: unable to locate plugin '%s'\n", plugin_type ENDFB(G);
    ok_raise(1);
  }

  // open file
  file_handle = plugin->open_file_read(fname, plugin_type, &natoms);

  if(!file_handle) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjectMolecule: plugin '%s' cannot open '%s'.\n", plugin_type, fname ENDFB(G);
    ok_raise(1);
  }

  // read atoms
  atoms = Calloc(molfile_atom_t, natoms);
  if (plugin->read_structure(file_handle, &optflags, atoms) != MOLFILE_SUCCESS) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjectMolecule: plugin '%s' failed to read atoms.\n", plugin_type ENDFB(G);
    ok_raise(1);
  }

  // Create ObjectMolecule
  ok_assert(1, I = ObjectMoleculeNew(G, false));
  I->Obj.Color = AtomInfoUpdateAutoColor(G);
  I->AtomInfo = VLACalloc(AtomInfoType, natoms);
  I->NAtom = natoms;

  // copy atom info
  for (int i = 0; i < natoms; i++) {
    AtomInfoType *ai = I->AtomInfo + i;
    molfile_atom_t *a = atoms + i;

    ai->rank = i;
    ai->id = i + 1;
    ai->b = a->bfactor;
    ai->q = a->occupancy;
    ai->vdw = a->radius;
    ai->partialCharge = a->charge;
    ai->alt[0] = a->altloc[0];

    strncpy(ai->segi, a->segid, cSegiLen);
    strncpy(ai->resn, a->resname, cResnLen);
    strncpy(ai->name, a->name, cAtomNameLen);
    if (a->atomicnumber > 0)
      atomicnumber2elem(ai->elem, a->atomicnumber);

    ai->chain = LexIdx(G, a->chain);
    ai->textType = LexIdx(G, a->type);

    ai->hetatm = 0;

    ai->resv = a->resid;
    snprintf(ai->resi, cResnLen, "%d%s", a->resid, a->insertion);

    ai->visRep = auto_show;

    AtomInfoAssignParameters(G, ai);
    AtomInfoAssignColors(G, ai);
  }

  // read coordinates
  while (/* true */ plugin->read_next_timestep != NULL) {
    ok_assert(1, cs = CoordSetNew(G));
    ok_assert(1, cs->Coord = VLAlloc(float, 3 * natoms));

    timestep.coords = cs->Coord;
    timestep.velocities = NULL;

    if (plugin->read_next_timestep(file_handle, natoms, &timestep) != MOLFILE_SUCCESS) {
      cs->fFree();
      break;
    }

    cs->Obj = I;
    cs->NIndex = natoms;
    cs->enumIndices();

    // append to object
    VLACheck(I->CSet, CoordSet*, I->NCSet);
    I->CSet[I->NCSet++] = cs;
  }

  // read bonds
  if (plugin->read_bonds &&
      plugin->read_bonds(file_handle, &nbonds, &from, &to, &order,
        &bondtype, &nbondtypes, &bondtypename) != MOLFILE_SUCCESS) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjectMolecule: plugin '%s' failed to read bonds.\n", plugin_type ENDFB(G);
    ok_raise(1);
  }

  // copy bonds
  if (nbonds) {
    I->NBond = nbonds;
    I->Bond = VLACalloc(BondType, nbonds);
    for (int i = 0; i < nbonds; i++) {
      BondTypeInit2(I->Bond + i, from[i] - 1, to[i] - 1,
          order ? (int) order[i] : 1);
    }
  } else if (I->NCSet) {
    ObjectMoleculeConnect(I, &I->NBond, &I->Bond, I->AtomInfo, I->CSet[0], true, -1);
  }

  // symmetry
  I->Symmetry = SymmetryNewFromTimestep(G, &timestep);

  // finalize
  ObjectMoleculeSort(I);
  ObjectMoleculeInvalidate(I, cRepAll, cRepInvAll, -1);
  ObjectMoleculeUpdateIDNumbers(I);
  ObjectMoleculeUpdateNonbonded(I);

  // merge
  if(origObj) {
    // TODO
  }

  SceneCountFrames(G);

ok_except1:
  // close
  if (plugin && file_handle)
    plugin->close_file_read(file_handle);

  if(atoms)
    free(atoms);

  return I;
}

/*
 * Load any object type with the given plugin. If obj_ptr's object type
 * doesn't match the plugin, the object will be deleted and a new one created
 * (not for trajectories).
 */
CObject * PlugIOManagerLoad(PyMOLGlobals * G, CObject ** obj_ptr,
    const char *fname, int state, int quiet, const char *plugin_type)
{
  CObject *obj = obj_ptr ? *obj_ptr : NULL;
  CPlugIOManager *manager = G->PlugIOManager;
  molfile_plugin_t *plugin;

  ok_assert(1, manager);
  plugin = find_plugin(manager, plugin_type);

  if (!plugin) {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " PlugIOManagerLoad: no plugin '%s'\n", plugin_type ENDFB(G);
    return NULL;
  }

  if (plugin->read_volumetric_data != NULL) {
    // maps

    if (obj && obj->type != cObjectMap) {
      ExecutiveDelete(G, obj->Name);
      obj = *obj_ptr = NULL;
    }

    return (CObject *) PlugIOManagerLoadVol(G, (ObjectMap *) obj,
        fname, state, quiet, plugin_type);

  } else if (plugin->read_structure != NULL) {
    // molecules

    if (obj
#if 0
        && obj->type != cObjectMolecule
#else
        // TODO no merge support yet, always delete existing object
#endif
        ) {
      ExecutiveDelete(G, obj->Name);
      obj = *obj_ptr = NULL;
    }

    return (CObject *) PlugIOManagerLoadMol(G, (ObjectMolecule *) obj,
        fname, state, quiet, plugin_type);

  } else if (plugin->read_next_timestep != NULL) {
    // trajectories

    float shift[] = {0.f, 0.f, 0.f};

    if (obj && obj->type != cObjectMolecule) {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " PlugIOManagerLoad: can't load trajectory into object '%s'\n", obj->Name ENDFB(G);
      return NULL;
    }

    PlugIOManagerLoadTraj(G, (ObjectMolecule *) obj,
        fname, state, 1, 1, 1, -1, -1, "all", 1, shift, quiet, plugin_type);
    return NULL;

  }

  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " PlugIOManagerLoad: '%s' doesn't provide any read function\n", plugin_type ENDFB(G);

ok_except1:
  return NULL;
}

#ifdef __cplusplus
}
#endif

#endif
