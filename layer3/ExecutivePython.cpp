#include "ExecutivePython.h"
#ifndef _PYMOL_NOPY
#include "ObjectAlignment.h"
#include "ObjectCGO.h"
#include "ObjectCallback.h"
#include "ObjectMap.h"
#include "P.h"

pymol::Result<> ExecutiveLoadObject(PyMOLGlobals* G,
    const char* oname, PyObject* model, int frame, int type, int finish,
    int discrete, int quiet, int zoom)
{
  ObjectNameType valid_name = "";
  pymol::CObject *origObj = NULL, *obj;
  OrthoLineType buf;
  buf[0] = 0;
  ExecutiveProcessObjectName(G, oname, valid_name);

  origObj = ExecutiveFindObjectByName(G, valid_name);

  /* TODO check for existing object of wrong type */

  switch (type) {
  case cLoadTypeChemPyModel: {
    if (origObj) {
      if (origObj->type != cObjectMolecule) {
        ExecutiveDelete(G, valid_name);
        origObj = NULL;
      } else {
        discrete = 1;
      }
    }
    PBlock(G); /*PBlockAndUnlockAPI(); */
    obj = ObjectMoleculeLoadChemPyModel(
        G, (ObjectMolecule*) origObj, model, frame, discrete);
    PUnblock(G); /*PLockAPIAndUnblock(); */
    if (!origObj) {
      if (obj) {
        ObjectSetName(obj, valid_name);
        ExecutiveManageObject(G, obj, zoom, quiet);
        if (frame < 0)
          frame = ((ObjectMolecule*) obj)->NCSet - 1;
        sprintf(buf,
            " CmdLoad: ChemPy-model loaded into object \"%s\", state %d.\n",
            valid_name, frame + 1);
      }
    } else if (origObj) {
      if (finish)
        ExecutiveUpdateObjectSelection(G, origObj);
      if (frame < 0)
        frame = ((ObjectMolecule*) origObj)->NCSet - 1;
      sprintf(buf,
          " CmdLoad: ChemPy-model appended into object \"%s\", state %d.\n",
          valid_name, frame + 1);
    }
    break;
  }
  case cLoadTypeChemPyBrick:
    if (origObj)
      if (origObj->type != cObjectMap) {
        ExecutiveDelete(G, valid_name);
        origObj = NULL;
      }
    PBlock(G); /*PBlockAndUnlockAPI(); */
    obj = ObjectMapLoadChemPyBrick(
        G, (ObjectMap*) origObj, model, frame, discrete, quiet);
    PUnblock(G); /*PLockAPIAndUnblock(); */
    if (!origObj) {
      if (obj) {
        ObjectSetName(obj, valid_name);
        ExecutiveManageObject(G, obj, zoom, quiet);
        sprintf(buf, " CmdLoad: chempy.brick loaded into object \"%s\"\n",
            valid_name);
      }
    } else if (origObj) {
      sprintf(buf, " CmdLoad: chempy.brick appended into object \"%s\"\n",
          valid_name);
    }
    break;
  case cLoadTypeChemPyMap:
    if (origObj)
      if (origObj->type != cObjectMap) {
        ExecutiveDelete(G, valid_name);
        origObj = NULL;
      }
    PBlock(G); /*PBlockAndUnlockAPI(); */
    obj = ObjectMapLoadChemPyMap(
        G, (ObjectMap*) origObj, model, frame, discrete, quiet);
    PUnblock(G); /*PLockAPIAndUnblock(); */
    if (!origObj) {
      if (obj) {
        ObjectSetName(obj, valid_name);
        ExecutiveManageObject(G, obj, zoom, quiet);
        sprintf(buf, " CmdLoad: chempy.map loaded into object \"%s\"\n",
            valid_name);
      }
    } else if (origObj) {
      sprintf(buf, " CmdLoad: chempy.map appended into object \"%s\"\n",
          valid_name);
    }
    break;
  case cLoadTypeCallback:
    if (origObj)
      if (origObj->type != cObjectCallback) {
        ExecutiveDelete(G, valid_name);
        origObj = NULL;
      }
    PBlock(G); /*PBlockAndUnlockAPI(); */
    obj = ObjectCallbackDefine(
        G, (ObjectCallback*) origObj, model, frame);
    PUnblock(G); /*PLockAPIAndUnblock(); */
    if (!origObj) {
      if (obj) {
        ObjectSetName(obj, valid_name);
        ExecutiveManageObject(G, obj, zoom, quiet);
        sprintf(buf, " CmdLoad: pymol.callback loaded into object \"%s\"\n",
            valid_name);
      }
    } else if (origObj) {
      sprintf(buf, " CmdLoad: pymol.callback appended into object \"%s\"\n",
          valid_name);
    }
    break;
  case cLoadTypeCGO:
    if (origObj)
      if (origObj->type != cObjectCGO) {
        ExecutiveDelete(G, valid_name);
        origObj = NULL;
      }
    PBlock(G); /*PBlockAndUnlockAPI(); */
    obj = ObjectCGODefine(G, (ObjectCGO*) origObj, model, frame);
    PUnblock(G); /*PLockAPIAndUnblock(); */
    if (!origObj) {
      if (obj) {
        ObjectSetName(obj, valid_name);
        ExecutiveManageObject(G, obj, zoom, quiet);
        sprintf(buf, " CmdLoad: CGO loaded into object \"%s\"\n", valid_name);
      }
    } else if (origObj) {
      sprintf(buf, " CmdLoad: CGO appended into object \"%s\"\n", valid_name);
    }
    break;
  }
  if (origObj && !quiet) {
    PRINTFB(G, FB_Executive, FB_Actions)
    "%s", buf ENDFB(G);
    OrthoRestorePrompt(G);
  }
  return {};
}

pymol::Result<> ExecutiveSetRawAlignment(PyMOLGlobals* G,
    pymol::zstring_view alnname, PyObject* raw, pymol::zstring_view guidename,
    int state, int quiet)
{

  ObjectMolecule* guide = nullptr;
  if (!guidename.empty()) {
    guide = ExecutiveFindObject<ObjectMolecule>(G, guidename.c_str());
  }

  if(!PyList_Check(raw)) {
    return pymol::make_error("alignment must be list");
  }

  auto n_cols = PyList_Size(raw);

  pymol::vla<int> align_vla(n_cols * 3);
  size_t vla_offset = 0;

  for(size_t c = 0; c < n_cols; ++c) {
    PyObject * col = PyList_GetItem(raw, c);

    if(!PyList_Check(col)) {
      return pymol::make_error("columns must be list");
    }

    auto n_idx = PyList_Size(col);

    for(size_t i = 0; i < n_idx; ++i) {
      const char * model;
      int index;

      PyObject * idx = PyList_GetItem(col, i);

      if(!PyArg_ParseTuple(idx, "si", &model, &index)) {
        return pymol::make_error("indices must be (str, int)");
      }

      auto mol = ExecutiveFindObject<ObjectMolecule>(G, model);

      if(!mol) {
        return pymol::make_error("object ", model, " not found");
      }

      if (!guide) {
        guide = mol;
      }

      if (index < 1 || mol->NAtom < index) {
        return pymol::make_error("index ('", model, ", ", index, ") out of range");
      }

      auto uid = AtomInfoCheckUniqueID(G, mol->AtomInfo + index - 1);
      *(align_vla.check(vla_offset++)) = uid;
    }

    *(align_vla.check(vla_offset++)) = 0;
  }

  align_vla.resize(vla_offset);

  // does alignment object already exist?
  auto cobj = ExecutiveFindObjectByName(G, alnname.c_str());
  if (cobj && cobj->type != cObjectAlignment) {
    ExecutiveDelete(G, cobj->Name);
    cobj = nullptr;
  }

  // create alignment object
  cobj = ObjectAlignmentDefine(G, (ObjectAlignment*) cobj,
      align_vla, state, true, guide, nullptr);

  // manage alignment object
  ObjectSetName(cobj, alnname.c_str());
  ExecutiveManageObject(G, cobj, 0, quiet);
  SceneInvalidate(G);

  // make available as selection FIXME find better solution
  cobj->update();
  return {};
}

pymol::Result<float> ExecutiveFitPairs(
    PyMOLGlobals* G, PyObject* list, int quiet)
{
  auto ln = PyObject_Length(list);
  if (!ln) {
    return pymol::make_error("No selections provided");
  }
  if (ln & 0x1) {
    return pymol::make_error(
        G, "FitPairs", "must supply an even number of selections.");
  }

  std::vector<SelectorTmp> word(ln);

  int a = 0;
  while (a < ln) {
    unique_PyObject_ptr item(PySequence_GetItem(list, a));
    auto tmp = SelectorTmp::make(G, PyString_AsString(item.get()));
    p_return_if_error(tmp);
    word[a] = std::move(tmp.result());
    a++;
  }
  return ExecutiveRMSPairs(G, word, 2, quiet);
}

#endif
