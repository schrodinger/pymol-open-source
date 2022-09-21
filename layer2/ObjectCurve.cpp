#include "ObjectCurve.h"
#include "CGO.h"
#include "Color.h"
#include "PConv.h"
#include "Picking.h"
#include "ShaderMgr.h"
#include "TTT.h"

ObjectCurveState::ObjectCurveState(PyMOLGlobals* G)
    : CObjectState(G)
{
}

void ObjectCurveState::addDefaultBezierSpline()
{
  if (splines.empty()) {
    splines.emplace_back();
    auto& firstSpline = splines.back();
    firstSpline.addBezierPoint();
  }
}

glm::vec3 ObjectCurveState::getPosition(float t) const
{
  const auto& spline = splines.front();
  return spline.getBezierPoint(t);
}

glm::vec3 ObjectCurveState::getNormalizedDirection(float t) const
{
  const auto& spline = splines.front();
  return spline.getFirstDerivative(t);
}

/**
 * @param cgo bezier cgo
 * @param curveState ObjectCurveState to validate
 * @param handleDiameter diameter of drawn validator spheres
 */

static void ValidateCGO(
    CGO& cgo, const ObjectCurveState& curveState, float handleDiameter)
{
  glm::vec3 redGLM = glm::vec3(1, 1, 1);
  auto red = glm::value_ptr(redGLM);
  int quality = 32;
  for (int i = 0; i < quality; i++) {
    float t = static_cast<float>(i) / static_cast<float>(quality);
    const auto pos = curveState.getPosition(t);
    cgo.add<cgo::draw::color>(red);
    cgo.add<cgo::draw::sphere>(glm::value_ptr(pos), handleDiameter);
  }
}

/**
 * Draws Handle for each control point
 * @param cgo bezier cgo
 * @param splineIndex index of spline to draw
 * @param pickIndex PickIndex
 * @param handleDiameter diameter of drawn handle spheres
 * @param ctrlPtPtr pointer to control point position
 * @param handlePointPtr pointer to handle point position
 */

static void DrawHandle(CGO& cgo, int splineIndex, int pickIndex,
    float handleDiameter, const float* ctrlPtPtr, const float* handlePtPtr)
{
  float handleLineRadius = 0.25f / 6.0f;
  glm::vec3 redGLM = glm::vec3(1, 0, 0);
  auto red = glm::value_ptr(redGLM);
  glm::vec3 whiteGLM = glm::vec3(1, 1, 1);
  auto white = glm::value_ptr(whiteGLM);
  CGOPickColor(&cgo, pickIndex, splineIndex);
  cgo.add<cgo::draw::color>(red);
  cgo.add<cgo::draw::sphere>(handlePtPtr, handleDiameter);

  cgo.add<cgo::draw::color>(white);
  cgo.add<cgo::draw::cylinder>(
      ctrlPtPtr, handlePtPtr, handleLineRadius, white, white);
}

void ObjectCurveState::updateRawCGO()
{
  rawCGO = nullptr;
  if (splines.empty()) {
    return;
  }
  auto& spline = splines.front();
  int splineIndex = 0; // TODO Iterate through splines
  const auto& bezierPoints = spline.getBezierPoints();
  rawCGO.reset(new CGO(G));
  float handleDiameter = 0.25f;
  for (int i = 1; i < bezierPoints.size(); ++i) {
    const auto& bezierA = bezierPoints[i - 1];
    const auto& bezierB = bezierPoints[i];
    rawCGO->add<cgo::draw::bezier>(glm::value_ptr(bezierA.control),
        glm::value_ptr(bezierA.rightHandle), glm::value_ptr(bezierB.leftHandle),
        glm::value_ptr(bezierB.control));
  }
  const glm::vec3 greenGLM(0, 1, 0.2);
  auto green = glm::value_ptr(greenGLM);

  for (int i = 0; i < bezierPoints.size(); ++i) {
    const auto& bezierPt = bezierPoints[i];
    auto pickIdx = i * 3;
    auto ctrlPt = glm::value_ptr(bezierPt.control);

    // draw control point
    CGOPickColor(rawCGO.get(), pickIdx, splineIndex);
    rawCGO->add<cgo::draw::color>(green);
    rawCGO->add<cgo::draw::sphere>(ctrlPt, handleDiameter);

    // draw Left handle -- Draw for all except first
    if (i != 0) {
      DrawHandle(*rawCGO, splineIndex, pickIdx + 1, handleDiameter, ctrlPt,
          glm::value_ptr(bezierPt.leftHandle));
    }

    // draw Right handle -- Draw for all except last
    if (i != bezierPoints.size() - 1) {
      DrawHandle(*rawCGO, splineIndex, pickIdx + 2, handleDiameter, ctrlPt,
          glm::value_ptr(bezierPt.rightHandle));
    }
  }

#if PYMOL_CURVE_VALIDATE
  ValidateCGO(*rawCGO, *this, handleDiameter);
#endif
}

static CGO* FilterCGO(PyMOLGlobals* G, const CGO* rawCGO)
{
  auto optCGO = pymol::make_unique<CGO>(G);
  CGO* allCylinders = nullptr;
  CGO* allBeziers = nullptr;
  CGO* allSpheres = nullptr;
  CGO* convertcgo = nullptr;
  if (CGOHasBezierOperations(rawCGO)) {
    CGO* allButBezier = new CGO(G);
    allBeziers = CGOOptimizeBezier(rawCGO);
    CGOFilterOutBezierOperationsInto(rawCGO, allButBezier);
    CGOStop(allButBezier);
    CGOFree(convertcgo);
    convertcgo = allButBezier;
  }
  if (CGOHasCylinderOperations(rawCGO)) {
    allCylinders = CGONew(G);
    CGOEnable(allCylinders, GL_CYLINDER_SHADER);
    CGO* newCGO =
        CGOConvertShaderCylindersToCylinderShader(rawCGO, allCylinders);
    allCylinders->free_append(newCGO);
    assert(newCGO == nullptr);
    CGODisable(allCylinders, GL_CYLINDER_SHADER);
    CGOStop(allCylinders);
    CGO* allButCylinders = CGONew(G);
    CGOFilterOutCylinderOperationsInto(rawCGO, allButCylinders);
    CGOStop(allButCylinders);
    CGOFree(convertcgo);
    convertcgo = allButCylinders;
  }
  if (CGOHasSphereOperations(rawCGO)) {
    CGO* allButSpheres = CGONew(G);
    allSpheres =
        CGOOptimizeSpheresToVBONonIndexed(rawCGO, 0, true, allButSpheres);
    if (allSpheres) {
      CGOFree(convertcgo);
      CGOStop(allButSpheres);
      convertcgo = allButSpheres;
    } else {
      CGOFree(allButSpheres);
    }
  }
  optCGO.reset(CGOSimplify(convertcgo));
  optCGO.reset(CGOOptimizeToVBONotIndexed(optCGO.get()));

  if (allBeziers) {
    optCGO->free_append(allBeziers);
  }
  if (allSpheres) {
    optCGO->free_append(allSpheres);
  }
  if (allCylinders) {
    optCGO->free_append(allCylinders);
  }
  return optCGO.release();
}

void ObjectCurveState::updateRenderCGO()
{
  if (renderCGO) { // We've already generated the renderCGO
    return;
  }
  if (!rawCGO) { // Raw CGO needs to be updated
    updateRawCGO();
    if (!rawCGO) {
      return;
    }
  }
  renderCGO.reset(FilterCGO(G, rawCGO.get()));
}

ObjectCurve::ObjectCurve(PyMOLGlobals* G)
    : pymol::CObject(G)
{
  type = cObjectCurve;

  m_states.emplace_back(G);
  auto& defaultState = m_states.back();
  defaultState.addDefaultBezierSpline();
}

void ObjectCurve::render(RenderInfo* info)
{
  ObjectPrepareContext(this, info);
  if (!(this->visRep & cRepCGOBit)) {
    return;
  }
  auto pass = info->pass;
  auto color = ColorGet(G, this->Color);

  if (info->ray) {
    return;
  }

  if (!G->HaveGUI || !G->ValidContext) {
    return;
  }

  for (auto state : StateIteratorV2(this, info->state)) {
    if (state >= m_states.size()) {
      continue;
    }
    auto& stateObj = m_states[state];
    if (info->pick) {
      PickContext context{};
      context.object = this;
      auto cgo = stateObj.renderCGO.get();
      CGORenderPicking(cgo, info, &context, Setting.get(), nullptr);
      continue;
    }
    if (pass != RenderPass::Antialias) {
      stateObj.updateRenderCGO();
      auto cgo = stateObj.renderCGO.get();
      if (cgo) {
        CGORender(cgo, color, this->Setting.get(), nullptr, info, nullptr);
      }
    }
  }
}

void ObjectCurve::update()
{
  for (auto& state : m_states) {
    state.renderCGO = nullptr;
  }
}

pymol::CObject* ObjectCurve::clone() const
{
  return new ObjectCurve(*this);
}

pymol::Result<pymol::BezierSplinePoint> ObjectCurve::getBezierPointByPick(
    const Picking& pick)
{
  assert(pick.context.state >= 0 && pick.context.state < m_states.size());
  const auto& state = m_states[pick.context.state];
  assert(pick.src.bond < state.splines.size());
  const auto& spline = state.splines[pick.src.bond];
  assert(pick.src.index < (spline.getBezierPoints().size() * 3));
  int bezPtIdx = pick.src.index / 3;
  return spline.getBezierPoints()[bezPtIdx];
}

pymol::Result<glm::vec3> ObjectCurve::getPositionByPick(const Picking& pick)
{
  auto bezPt = getBezierPointByPick(pick);
  if (!bezPt) {
    return bezPt.error();
  }
  int handleIdx = pick.src.index % 3;
  glm::vec3 pt;
  switch (handleIdx) {
  case 0:
    pt = bezPt->control;
    break;
  case 1:
    pt = bezPt->leftHandle;
    break;
  case 2:
    pt = bezPt->rightHandle;
    break;
  default:
    return pymol::make_error("Invalid handle index");
    break;
  }
  if (TTTFlag) {
    pt = pymol::TTT::from_pymol_2_legacy(TTT).transform(pt);
  }
  return pt;
}

pymol::Result<> ObjectCurve::setPositionByPick(
    const Picking& pick, const glm::vec3& newPos)
{
  assert(pick.context.state >= 0 && pick.context.state < m_states.size());
  auto& state = m_states[pick.context.state];
  assert(pick.src.bond < state.splines.size());
  auto& spline = state.splines[pick.src.bond];
  assert(pick.src.index < (spline.getBezierPoints().size() * 3));

  int bezPtIdx = pick.src.index / 3;
  int handleIdx = pick.src.index % 3;
  auto& bezPt = spline.getBezierPoints()[bezPtIdx];
  switch (handleIdx) {
  case 0: {
    auto delta = newPos - bezPt.control;
    bezPt.control += delta;
    bezPt.leftHandle += delta;
    bezPt.rightHandle += delta;
  } break;
  case 1: {
    bezPt.leftHandle = newPos;
    auto handleVec = bezPt.leftHandle - bezPt.control;
    if (bezPt.mode == pymol::BezierControlPointMode::ALIGNED) {
      bezPt.rightHandle = bezPt.control - handleVec;
    }
  } break;
  case 2: {
    bezPt.rightHandle = newPos;
    auto handleVec = bezPt.rightHandle - bezPt.control;
    if (bezPt.mode == pymol::BezierControlPointMode::ALIGNED) {
      bezPt.leftHandle = bezPt.control - handleVec;
    }
  } break;
  default:
    break;
  }
  state.rawCGO = nullptr;
  state.renderCGO = nullptr;
  return {};
}

pymol::Result<pymol::BezierSpline*> ObjectCurve::getBezierSplineByPick(
    const Picking& pick)
{
  assert(pick.context.state >= 0 && pick.context.state < m_states.size());
  auto& state = m_states[pick.context.state];
  assert(pick.src.bond < state.splines.size());
  auto& spline = state.splines[pick.src.bond];
  return &spline;
}

glm::vec3 ObjectCurve::getPosition(float t) const
{
  int currentState = 0;
  const auto& state = m_states[currentState];
  auto rawPos = state.getPosition(t);
  if (!TTTFlag) {
    return rawPos;
  }
  auto ttt = pymol::TTT::from_pymol_2_legacy(TTT);
  return ttt.transform(rawPos);
}

glm::vec3 ObjectCurve::getNormalizedDirection(float t) const
{
  int currentState = 0;
  const auto& state = m_states[currentState];
  return state.getNormalizedDirection(t);
}

void ObjectCurve::invalidate(cRep_t rep, cRepInv_t level, int state)
{
  for (auto& state : m_states) {
    state.rawCGO = nullptr;
    state.renderCGO = nullptr;
  }
}

////////////////////////////
/////////SERIALIZE//////////
////////////////////////////

static PyObject* BezierPointAsPyList(const pymol::BezierSplinePoint& pt)
{
  auto numFloatsPerPoint = 10; // 3 vec3s + 1 enum
  auto bezierList = PyList_New(numFloatsPerPoint);
  PyList_SetItem(bezierList, 0, PyFloat_FromDouble(pt.control[0]));
  PyList_SetItem(bezierList, 1, PyFloat_FromDouble(pt.control[1]));
  PyList_SetItem(bezierList, 2, PyFloat_FromDouble(pt.control[2]));
  PyList_SetItem(bezierList, 3, PyFloat_FromDouble(pt.leftHandle[0]));
  PyList_SetItem(bezierList, 4, PyFloat_FromDouble(pt.leftHandle[1]));
  PyList_SetItem(bezierList, 5, PyFloat_FromDouble(pt.leftHandle[2]));
  PyList_SetItem(bezierList, 6, PyFloat_FromDouble(pt.rightHandle[0]));
  PyList_SetItem(bezierList, 7, PyFloat_FromDouble(pt.rightHandle[1]));
  PyList_SetItem(bezierList, 8, PyFloat_FromDouble(pt.rightHandle[2]));
  PyList_SetItem(bezierList, 9, PyInt_FromLong(static_cast<int>(pt.mode)));
  return PConvAutoNone(bezierList);
}

static PyObject* BezierSplineAsPyList(const pymol::BezierSpline& spline)
{
  const auto& bezPts = spline.getBezierPoints();
  auto result = PyList_New(bezPts.size());
  for (int i = 0; i < bezPts.size(); i++) {
    PyList_SetItem(result, i, BezierPointAsPyList(bezPts[i]));
  }
  return PConvAutoNone(result);
}

PyObject* ObjectCurveState::asPyList() const
{
  auto result = PyList_New(splines.size());
  for (int i = 0; i < splines.size(); i++) {
    PyList_SetItem(result, i, BezierSplineAsPyList(splines[i]));
  }

  return PConvAutoNone(result);
}

PyObject* ObjectCurve::statesAsPyList() const
{
  auto result = PyList_New(m_states.size());
  for (int a = 0; a < m_states.size(); a++) {
    PyList_SetItem(result, a, m_states[a].asPyList());
  }
  return PConvAutoNone(result);
}

PyObject* ObjectCurve::asPyList() const
{
  auto result = PyList_New(2);
  PyList_SetItem(result, 0, ObjectAsPyList(this));
  PyList_SetItem(result, 1, statesAsPyList());
  return PConvAutoNone(result);
}

////////////////////////////
////////DESERIALIZE/////////
////////////////////////////

pymol::Result<pymol::BezierSplinePoint> BezierSplineFromPyList(
    PyObject* serializedList)
{
  pymol::BezierSplinePoint pt;
  if (!PyList_Check(serializedList)) {
    return pymol::make_error("BezierSplinePoint: Not a list");
  }
  assert(PyList_Size(serializedList) == 10); // 3 vec3s + 1 enum
  auto to_float = [&](int idx) {
    return static_cast<float>(
        PyFloat_AsDouble(PyList_GetItem(serializedList, idx)));
  };
  auto to_glm_vec3 = [&](int idx) {
    return glm::vec3(to_float(idx), to_float(idx + 1), to_float(idx + 2));
  };
  auto to_handle_mode = [&](int idx) {
    return static_cast<pymol::BezierControlPointMode>(
        PyInt_AsLong(PyList_GetItem(serializedList, idx)));
  };
  pt.control = to_glm_vec3(0);
  pt.leftHandle = to_glm_vec3(3);
  pt.rightHandle = to_glm_vec3(6);
  pt.mode = to_handle_mode(9);
  return pt;
}

ObjectCurveState::ObjectCurveState(PyMOLGlobals* G, PyObject* serializedList)
    : CObjectState(G)
{
  if (!PyList_Check(serializedList)) {
    printf("ObjectCurveState: Could not deserialize list\n");
    return;
  }
  int newSplinesSize = PyList_Size(serializedList);
  for (int a = 0; a < newSplinesSize; a++) {
    auto splinePyList = CPythonVal_PyList_GetItem(I->G, serializedList, a);
    splines.emplace_back();
    auto& newSpline = splines.back();
    int newPointsListSize = PyList_Size(splinePyList);
    for (int bp = 0; bp < newPointsListSize; bp++) {
      auto newPointList = CPythonVal_PyList_GetItem(I->G, splinePyList, bp);
      if (auto splineListResult = BezierSplineFromPyList(newPointList)) {
        newSpline.addBezierPoint(*splineListResult);
      }
    }
  }
}

pymol::Result<> ObjectCurve::statesFromPyList(PyObject* serializedList)
{
  if (!PyList_Check(serializedList)) {
    return pymol::make_error("Curve States: Invalid PyList");
  }
  int newStatesSize = PyList_Size(serializedList);
  for (int a = 0; a < newStatesSize; a++) {
    auto statePyList = CPythonVal_PyList_GetItem(I->G, serializedList, a);
    m_states.emplace_back(G, statePyList);
    CPythonVal_Free(statePyList);
  }
  return {};
}

ObjectCurve::ObjectCurve(PyMOLGlobals* G, PyObject* serializedList)
    : pymol::CObject(G)
{
  bool ok = true;
  if (ok) {
    auto val = CPythonVal_PyList_GetItem(G, serializedList, 0);
    ok = ObjectFromPyList(G, val, this);
    CPythonVal_Free(val);
  }
  if (ok) {
    auto val = CPythonVal_PyList_GetItem(G, serializedList, 1);
    statesFromPyList(val);
    CPythonVal_Free(val);
  }
}
