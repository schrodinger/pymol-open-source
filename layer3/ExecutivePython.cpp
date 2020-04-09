#include "ExecutivePython.h"
#include "ObjectCGO.h"
#include "ObjectCallback.h"
#include "ObjectMap.h"
#include "P.h"

pymol::Result<> ExecutiveLoadObject(PyMOLGlobals* G,
    const char* oname, PyObject* model, int frame, int type, int finish,
    int discrete, int quiet, int zoom)
{
  ObjectNameType valid_name = "";
  CObject *origObj = NULL, *obj;
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
#ifndef _PYMOL_NO_UNDO
#endif
    obj = (CObject*) ObjectMoleculeLoadChemPyModel(
        G, (ObjectMolecule*) origObj, model, frame, discrete);
    PUnblock(G); /*PLockAPIAndUnblock(); */
    if (!origObj) {
      if (obj) {
        ObjectSetName(obj, valid_name);
#ifndef _PYMOL_NO_UNDO
#endif
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
    obj = (CObject*) ObjectMapLoadChemPyBrick(
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
    obj = (CObject*) ObjectMapLoadChemPyMap(
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
    obj = (CObject*) ObjectCallbackDefine(
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
    obj = (CObject*) ObjectCGODefine(G, (ObjectCGO*) origObj, model, frame);
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
