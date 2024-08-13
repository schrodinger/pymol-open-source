#include "Bezier.h"
#include "PyMOLObject.h"
#include "Result.h"

#include <glm/vec3.hpp>
#include <vector>

struct ObjectCurveState : public CObjectState {
  ObjectCurveState(PyMOLGlobals* G);

  /**
   * Deserializes pylists into ObjectCurveState
   * @param serializedList serialized ObjectCurveState
   */
  ObjectCurveState(PyMOLGlobals* G, PyObject* serializedList);

  /**
   * Adds spline to object
   */
  void addDefaultBezierSpline();

  /**
   * @param t parameter t across Bezier Spline
   * @return interpolated position along curve object
   */
  glm::vec3 getPosition(float t) const;

  /**
   * @param t parameter t across Bezier Spline
   * @return direction at t
   */
  glm::vec3 getNormalizedDirection(float t) const;

  /**
   * updates the render CGO
   */
  void updateRenderCGO();

  /**
   * updates the raw bezier CGO
   */
  void updateRawCGO();

  /**
   * Serializes ObjectCurveState into python lists
   */
  PyObject* asPyList() const;

  std::vector<pymol::BezierSpline> splines;
  pymol::cache_ptr<CGO> rawCGO;
  pymol::cache_ptr<CGO> renderCGO;
};

class ObjectCurve : public pymol::CObject
{
public:
  ObjectCurve(PyMOLGlobals* G);

  /**
   * Deserializes pylists into ObjectCurve
   * @param serializedList serialized ObjectCurve
   */
  ObjectCurve(PyMOLGlobals* G, PyObject* serializedList);

  /**
   * Creates copy of this curve object
   */
  pymol::CObject* clone() const override;

  /**
   * Renders this object
   * @info render information context
   */
  void render(RenderInfo* info) override;

  /**
   * Updates the object
   */
  void update() override;

  /**
   * Invalidates parts of the object
   * @param rep representation[unused]
   * @param level invalidation level [unused]
   * @param state state to invalidate
   */
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;

  /**
   * For picking, we hijack:
   * Picking::src::index => bezier point index in spline
   * - 0: Control Point
   * - 1: Left Handle
   * - 2: Right Handle
   * Picking::src::bond => spline index
   */

  /**
   * @param pick Picking information (maps picking to object details)
   * @return position from picking information
   */
  pymol::Result<glm::vec3> getPositionByPick(const Picking& pick);

  /**
   * @param pick Picking information (maps picking to object details)
   * @param newPos new position of picked bezier point
   */
  pymol::Result<> setPositionByPick(
      const Picking& pick, const glm::vec3& newPos);

  /**
   * @param pick Picking information (maps picking to object details)
   * @return pointer to Bezier spline based on picking info
   */
  pymol::Result<pymol::BezierSpline*> getBezierSplineByPick(
      const Picking& pick);

  /**
   * Serializes ObjectCurve into python lists
   */
  PyObject* asPyList() const;

  /**
   * Serializes managed ObjectCurveStates into python lists
   */
  PyObject* statesAsPyList() const;

  /**
   * Deserializes ObjectCurveStates from python lists
   */
  pymol::Result<> statesFromPyList(PyObject* serializedList);

  /**
   * @param t parameter t across Bezier Spline
   * @return interpolated position along curve object
   */
  glm::vec3 getPosition(float t) const;

  /**
   * @param t parameter t across Bezier Spline
   * @return direction at t
   */
  glm::vec3 getNormalizedDirection(float t) const;

private:
  std::vector<ObjectCurveState> m_states;

  /**
   * @param pick Picking information (maps picking to object details)
   * @return BezierSplinePoint based on pick
   */
  pymol::Result<pymol::BezierSplinePoint> getBezierPointByPick(
      const Picking& pick);
};
