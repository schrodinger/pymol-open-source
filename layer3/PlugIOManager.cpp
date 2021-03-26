
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

#include <vector>

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
#include "Lex.h"
#include "CGO.h"
#include "ObjectCGO.h"
#include "Util.h"

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
                          int stop, int max, const char *sele, int image,
                          const float *shift, int quiet, const char *plugin_type)
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

pymol::CObject * PlugIOManagerLoad(PyMOLGlobals * G, pymol::CObject ** obj_ptr,
    const char *fname, int state, int quiet, const char *plugin_type, int mask)
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
  if((I = (G->PlugIOManager = pymol::calloc<CPlugIOManager>(1)))) {
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

static CSymmetry* SymmetryNewFromTimestep(
    PyMOLGlobals* G, molfile_timestep_t* ts);

int PlugIOManagerLoadTraj(PyMOLGlobals * G, ObjectMolecule * obj,
                          const char *fname, int frame,
                          int interval, int average, int start,
                          int stop, int max, const char *sele, int image,
                          const float *shift, int quiet, const char *plugin_type)
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

  if (obj->DiscreteFlag) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " %s: Discrete objects not supported\n", __func__ ENDFB(G);
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
      } else if(natoms != obj->NAtom || (cs && cs->NIndex != natoms)) {
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
        ok_assert(1, cs->Coord = pymol::vla<float>(3 * natoms));

        cs->Obj = obj;
        cs->NIndex = natoms;
        cs->enumIndices();
      }

      auto xref = LoadTrajSeleHelper(obj, cs, sele);

      auto coordbuf = std::vector<float>(natoms * 3);
      timestep.coords = coordbuf.data();

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
                    // TODO this doesn't make any sense
                    float* fp = timestep.coords;
                    for (int i = 0; i < natoms; ++i) {
                      *(fp++) /= n_avg;
                      *(fp++) /= n_avg;
                      *(fp++) /= n_avg;
                    }
                  }
                  /* add new coord set */

                  for (int i = 0; i < natoms; ++i) {
                    int idx = xref ? xref[i] : i;
                    if (idx >= 0) {
                      assert(idx < cs->NIndex);
                      copy3(timestep.coords + 3 * i, cs->coordPtr(idx));
                    }
                  }

                  cs->invalidateRep(cRepAll, cRepInvRep);
                  if(frame < 0) frame = obj->NCSet;
                  if(!obj->NCSet) zoom_flag = true;

		  /* make sure we have room for 'frame' CoordSet*'s in obj->CSet */
		  /* TODO: TEST this function */
                  VLACheck(obj->CSet, CoordSet*, frame); /* was CoordSet* */
		  /* bump the object's state count */
                  if(obj->NCSet <= frame) obj->NCSet = frame + 1;
		  /* if there's data in this state's coordset, emtpy it */
                  delete obj->CSet[frame];
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

                  // symmetry
                  cs->Symmetry.reset(SymmetryNewFromTimestep(G, &timestep));

                  if((stop > 0 && cnt >= stop) || (max > 0 && ncnt >= max)) {
                    cs = NULL;
                    break;
                  }

                  frame++;
                  /* make a new cs */
                  cs = CoordSetCopy(cs);        /* otherwise, we need a place to put the next set */
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
        delete cs;
        SceneChanged(G);
        SceneCountFrames(G);
        if(zoom_flag)
          if(SettingGetGlobal_i(G, cSetting_auto_zoom)) {
            ExecutiveWindowZoom(G, obj->Name, 0.0, -1, 0, 0, quiet);        /* auto zoom (all states) */
          }

        auto const defer_limit = SettingGet<int>(G, cSetting_auto_defer_builds);
        if (defer_limit >= 0                   //
            && obj->getNFrame() >= defer_limit //
            && SettingGet<int>(G, cSetting_defer_builds_mode) <= 0) {
          PRINTFB(G, FB_ObjectMolecule, FB_Details)
          " ObjectMolecule-Details: Enabling defer_builds_mode\n" ENDFB(G);
          SettingSet(G, cSetting_defer_builds_mode, 3);
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

          ok_assert(1, datablock = pymol::malloc<float>(size));

          if(plugin->read_volumetric_data(file_handle, i, datablock, NULL) != MOLFILE_SUCCESS) {
            PRINTFB(G, FB_ObjectMolecule, FB_Errors)
              " PlugIOManager: read_volumetric_data failed\n" ENDFB(G);
            ok_raise(1);
          }

          {
            ObjectMapState *ms = NULL;

            if(!obj)
              ok_assert(1, obj = new ObjectMap(G));

            if(state < 0)
              state = obj->State.size();
            if(obj->State.size() <= state) {
              VecCheckEmplace (obj->State, state, G);
            }
            ms = &obj->State[state];

            ms->FDim[0] = v->xsize;
            ms->FDim[1] = v->ysize;
            ms->FDim[2] = v->zsize;
            ms->FDim[3] = 3;

            ms->Grid = std::vector<float>(3);
            ms->Dim = std::vector<int>(3);
            ms->Origin = std::vector<float>(3, 0.0f);
            ms->Range = std::vector<float>(3);

            float axes33f[9];
            float originf[3];
            copy3(v->xaxis, axes33f + 0);
            copy3(v->yaxis, axes33f + 3);
            copy3(v->zaxis, axes33f + 6);
            copy3(v->origin, originf);

            // check if inverted volume (TetsurfVolume would get normals wrong)
            bool inverted = determinant33f(axes33f) < 0;
            if (inverted) {
              // flip the z-axis
              add3f(axes33f + 6, originf, originf);
              invert3f(axes33f + 6);
            }

            // special case: orthogonal & cartesian-aligned
            // -> don't use State.Matrix which causes trouble for e.g. CCP4 export
            // (non-standard header) and ObjectSlice (incomplete implementation)
            bool aligned_axes =
              axes33f[0] > 0 && axes33f[4] > 0 && axes33f[8] > 0 &&
              fabs(axes33f[1]) <= R_SMALL4 && fabs(axes33f[2]) <= R_SMALL4 &&
              fabs(axes33f[3]) <= R_SMALL4 && fabs(axes33f[5]) <= R_SMALL4 &&
              fabs(axes33f[6]) <= R_SMALL4 && fabs(axes33f[7]) <= R_SMALL4;

            // set corners to a unit cube, and manage world space with the state matrix
            if (!aligned_axes) {
              PRINTFB(G, FB_ObjectMap, FB_Warnings)
                " Warning: Axes not aligned, using map-state matrix\n" ENDFB(G);

              double m44d[16];

              if(ms->Matrix.empty())
                ms->Matrix = std::vector<double>(16);

              // state matrix transformation
              copy33f44d(axes33f, m44d);
              copy3(originf, m44d + 12);
              transpose44d44d(m44d, ms->Matrix.data());
            }

            // axis and corner stuff in a unit cube
            {
              // prime min+max
              zero3f(ms->ExtentMin);
              ones3f(ms->ExtentMax);
              ones3f(ms->Range.data());

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

            if (aligned_axes) {
              ms->Grid[0] = axes33f[0] / (ms->FDim[0] - 1);
              ms->Grid[1] = axes33f[4] / (ms->FDim[1] - 1);
              ms->Grid[2] = axes33f[8] / (ms->FDim[2] - 1);

              for(int a = 0; a < 3; a++) {
                ms->Origin[a] = originf[a];
                ms->Range[a] = ms->Grid[a] * (ms->Dim[a] - 1);
                ms->ExtentMin[a] = ms->Origin[a];
                ms->ExtentMax[a] = ms->Origin[a] + ms->Grid[a] * ms->Max[a];
              }

              // corners
              for(int c = 0, d = 0; c < ms->FDim[2]; c += ms->FDim[2] - 1) {
                float v[3];
                v[2] = ms->Origin[2] + ms->Grid[2] * c;
                for(int b = 0; b < ms->FDim[1]; b += ms->FDim[1] - 1) {
                  v[1] = ms->Origin[1] + ms->Grid[1] * b;
                  for(int a = 0; a < ms->FDim[0]; a += ms->FDim[0] - 1, ++d) {
                    v[0] = ms->Origin[0] + ms->Grid[0] * a;
                    copy3f(v, ms->Corner + 3 * d);
                  }
                }
              }
            }

            // field
            ms->Field.reset(new Isofield(G, ms->FDim));
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
                int cc = inverted ? (ms->FDim[2] - c - 1) : c;
                for(b = 0; b < ms->FDim[1]; b++) {
                  for(a = 0; a < ms->FDim[0]; a++) {
                    F3(ms->Field->data, a, b, cc) = *(data_ptr++);
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

static CSymmetry * SymmetryNewFromTimestep(PyMOLGlobals * G, molfile_timestep_t * ts)
{
  CSymmetry * symm = NULL;
  ok_assert(1,
      ts->A > 0.f && ts->B > 0.f && ts->C > 0.f &&
      ts->alpha > 0.f && ts->beta > 0.f && ts->gamma > 0.f);
  ok_assert(1, symm = new CSymmetry(G));
  symm->Crystal.setDims(ts->A, ts->B, ts->C);
  symm->Crystal.setAngles(ts->alpha, ts->beta, ts->gamma);
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
  auto literal_names = SettingGet<bool>(G, cSetting_pdb_literal_names);

  memset(&timestep, 0, sizeof(molfile_timestep_t));

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
  atoms = pymol::calloc<molfile_atom_t>(natoms);
  if (plugin->read_structure(file_handle, &optflags, atoms) != MOLFILE_SUCCESS) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " ObjectMolecule: plugin '%s' failed to read atoms.\n", plugin_type ENDFB(G);
    ok_raise(1);
  }

  // Create ObjectMolecule
  ok_assert(1, I = new ObjectMolecule(G, false));
  I->Color = AtomInfoUpdateAutoColor(G);
  VLASize(I->AtomInfo, AtomInfoType, natoms);
  I->NAtom = natoms;

  // copy atom info
  for (int i = 0; i < natoms; i++) {
    AtomInfoType *ai = I->AtomInfo + i;
    molfile_atom_t *a = atoms + i;

    if (!literal_names) {
      UtilCleanStr(a->segid);
      UtilCleanStr(a->chain);
      UtilCleanStr(a->resname);
      UtilCleanStr(a->name);
    }

    ai->rank = i;
    ai->id = i + 1;
    ai->b = a->bfactor;
    ai->q = a->occupancy;
    ai->vdw = a->radius;
    ai->partialCharge = a->charge;
    ai->alt[0] = a->altloc[0];

    ai->segi = LexIdx(G, a->segid);
    ai->resn = LexIdx(G, a->resname);
    ai->name = LexIdx(G, a->name);
    if (a->atomicnumber > 0)
      atomicnumber2elem(ai->elem, a->atomicnumber);

    ai->chain = LexIdx(G, a->chain);
    ai->textType = LexIdx(G, a->type);

    ai->hetatm = 0;

    ai->resv = a->resid;
    ai->setInscode(a->insertion[0]);

    ai->visRep = auto_show;

    AtomInfoAssignParameters(G, ai);
    AtomInfoAssignColors(G, ai);
  }

  // read coordinates
  while (/* true */ plugin->read_next_timestep != NULL) {
    ok_assert(1, cs = CoordSetNew(G));
    ok_assert(1, cs->Coord = pymol::vla<float>(3 * natoms));

    timestep.coords = cs->Coord.data();
    timestep.velocities = NULL;

    if (plugin->read_next_timestep(file_handle, natoms, &timestep) != MOLFILE_SUCCESS) {
      delete cs;
      break;
    }

    cs->Obj = I;
    cs->NIndex = natoms;
    cs->enumIndices();

    // symmetry
    cs->Symmetry.reset(SymmetryNewFromTimestep(G, &timestep));

    // append to object
    VLACheck(I->CSet, CoordSet*, I->NCSet);
    I->CSet[I->NCSet++] = cs;
  }

  // topology-only (with template coord set)
  if (!I->NCSet) {
    ok_assert(1, cs = CoordSetNew(G));
    ok_assert(1, cs->Coord = pymol::vla<float>(3 * natoms));

    cs->Obj = I;
    cs->NIndex = natoms;
    cs->enumIndices();

    I->CSTmpl = cs;
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
    I->Bond = pymol::vla<BondType>(nbonds);
    for (int i = 0; i < nbonds; i++) {
      BondTypeInit2(I->Bond + i, from[i] - 1, to[i] - 1,
          order ? (int) order[i] : 1);
    }
  } else if (I->NCSet) {
    ObjectMoleculeConnect(I, I->CSet[0]);
  }

  // finalize
  I->invalidate(cRepAll, cRepInvAll, -1);
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
    mfree(atoms);

  return I;
}

static void cgo_check_beginend(int type, int & currenttype, CGO *& cgo) {
  if (currenttype != type) {
    if (currenttype)
      CGOEnd(cgo);
    if (type)
      CGOBegin(cgo, type);
    currenttype = type;
  }
}

static
ObjectCGO *PlugIOManagerLoadGraphics(PyMOLGlobals * G, ObjectCGO *origObj,
    const char *fname, int state, int quiet, const char *plugin_type)
{
  CPlugIOManager *manager = G->PlugIOManager;
  void *file_handle = NULL;
  molfile_plugin_t * plugin = NULL;
  const molfile_graphics_t * graphics = NULL;
  int nelem = 0;
  int beginend = 0;
  CGO *cgo = NULL;
  ObjectCGO *I = NULL;

  ok_assert(1, manager);
  plugin = find_plugin(manager, plugin_type);

  if (!plugin) {
    PRINTFB(G, FB_ObjectCGO, FB_Errors)
      " ObjectCGO: unable to locate plugin '%s'\n", plugin_type ENDFB(G);
    ok_raise(1);
  }

  // open file
  file_handle = plugin->open_file_read(fname, plugin_type, &nelem /* dummy */);

  if(!file_handle) {
    PRINTFB(G, FB_ObjectCGO, FB_Errors)
      " ObjectCGO: plugin '%s' cannot open '%s'.\n", plugin_type, fname ENDFB(G);
    ok_raise(1);
  }

  // read geometry
  if (plugin->read_rawgraphics(file_handle, &nelem, &graphics) != MOLFILE_SUCCESS) {
    PRINTFB(G, FB_ObjectCGO, FB_Errors)
      " ObjectCGO: plugin '%s' failed to read graphics.\n", plugin_type ENDFB(G);
    ok_raise(1);
  }

  cgo = CGONew(G);

#define CHECK_BEGINEND(type) cgo_check_beginend(type, beginend, cgo)

  // translate to CGO
  for (auto g = graphics, g_end = graphics + nelem; g != g_end; ++g) {
    auto g_current = g;
    const float * tnormals = NULL;
    const float * tcolors = NULL;

    switch (g->type) {
      case MOLFILE_POINT:
        break;
      case MOLFILE_TRINORM:
        // followed by normals
      case MOLFILE_TRICOLOR:
        // followed by normals and color
        if (g + 1 != g_end && g[1].type == MOLFILE_NORMS) {
          tnormals = (++g)->data;
        }
        if (g_current->type == MOLFILE_TRICOLOR &&
            g + 1 != g_end && g[1].type == MOLFILE_COLOR) {
          tcolors = (++g)->data;
        }
      case MOLFILE_TRIANGLE:
        CHECK_BEGINEND(GL_TRIANGLES);
        for (int i = 0; i < 9; i += 3) {
          if (tnormals)
            CGONormalv(cgo, tnormals + i);
          if (tcolors)
            CGOColorv(cgo, tcolors + i);
          CGOVertexv(cgo, g_current->data + i);
        }
        break;
      case MOLFILE_NORMS:
        CGONormalv(cgo, g->data);
        break;
      case MOLFILE_LINE:
        CHECK_BEGINEND(GL_LINES);
        CGOVertexv(cgo, g->data + 0);
        CGOVertexv(cgo, g->data + 3);
        break;
      case MOLFILE_CYLINDER:
        CHECK_BEGINEND(0);
        {
          float axis[3];
          subtract3f(g->data + 3, g->data, axis);
          CGOShaderCylinder(cgo, g->data, axis, g->size, 0);
        }
        break;
      case MOLFILE_CAPCYL:
        break;
      case MOLFILE_CONE:
        break;
      case MOLFILE_SPHERE:
        CHECK_BEGINEND(0);
        CGOSphere(cgo, g->data, g->size);
        break;
      case MOLFILE_TEXT:
        break;
      case MOLFILE_COLOR:
        CGOColorv(cgo, g->data);
        break;
    }
  }

  CHECK_BEGINEND(0);
  CGOStop(cgo);

  // Create ObjectCGO
  ok_assert(1, I = ObjectCGOFromCGO(G, NULL, cgo, state));

  // default is cgo_lighting=0 when loading CGOs without normals
  SettingSet(cSetting_cgo_lighting, 1, (pymol::CObject *)I);

ok_except1:
  // close
  if (plugin && file_handle)
    plugin->close_file_read(file_handle);

  if (!I)
    CGOFree(cgo);

  return I;
}

/**
 * Load any object type with the given plugin. If obj_ptr's object type
 * doesn't match the plugin, the object will be deleted and a new one created
 * (not for trajectories).
 */
pymol::CObject * PlugIOManagerLoad(PyMOLGlobals * G, pymol::CObject ** obj_ptr,
    const char *fname, int state, int quiet, const char *plugin_type, int mask)
{
  pymol::CObject *obj = obj_ptr ? *obj_ptr : NULL;
  CPlugIOManager *manager = G->PlugIOManager;
  molfile_plugin_t *plugin;

  ok_assert(1, manager);
  plugin = find_plugin(manager, plugin_type);

  if (!plugin) {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " PlugIOManagerLoad: no plugin '%s'\n", plugin_type ENDFB(G);
    return NULL;
  }

  if (!mask)
    mask = cPlugIOManager_any;

  if ((mask & cPlugIOManager_vol) && plugin->read_volumetric_data) {
    // maps

    if (obj && obj->type != cObjectMap) {
      ExecutiveDelete(G, obj->Name);
      obj = *obj_ptr = NULL;
    }

    return PlugIOManagerLoadVol(G, (ObjectMap *) obj,
        fname, state, quiet, plugin_type);

  } else if ((mask & cPlugIOManager_mol) && plugin->read_structure) {
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

    return PlugIOManagerLoadMol(G, (ObjectMolecule *) obj,
        fname, state, quiet, plugin_type);

  } else if ((mask & cPlugIOManager_traj) && plugin->read_next_timestep) {
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

  } else if ((mask & cPlugIOManager_graphics) && plugin->read_rawgraphics) {
    // geometry (CGO)

    if (obj) {
      // no merge support
      ExecutiveDelete(G, obj->Name);
      obj = *obj_ptr = NULL;
    }

    return PlugIOManagerLoadGraphics(G, (ObjectCGO *) obj,
        fname, state, quiet, plugin_type);
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

/**
 * Find a plugin by filename extension
 *
 * @param ext File extension
 * @param mask plugin needs to read any content (0), structure (1), trajectory (2) or map (4)
 */
const char * PlugIOManagerFindPluginByExt(PyMOLGlobals * G, const char * ext, int mask) {
#ifdef _PYMOL_VMD_PLUGINS
  CPlugIOManager *I = G->PlugIOManager;

  if (!mask)
    mask = cPlugIOManager_any;

  for (auto it = I->PluginVLA, it_end = it + I->NPlugin; it != it_end; ++it) {
    const molfile_plugin_t * p = *it;

    if (WordMatchCommaExact(G, p->filename_extension, ext, true) >= 0)
      continue;

    if (((mask & cPlugIOManager_mol) && p->read_structure) ||
        ((mask & cPlugIOManager_traj) && p->read_next_timestep) ||
        ((mask & cPlugIOManager_graphics) && p->read_rawgraphics) ||
        ((mask & cPlugIOManager_vol) && p->read_volumetric_data))
      return p->name;
  }
#endif

  return NULL;
}
