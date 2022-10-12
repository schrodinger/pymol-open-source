#include "Bezier.h"

#include "pymol/algorithm.h"

#include <glm/gtx/io.hpp>

namespace pymol
{
//////////////////////////
///////PUBLIC API/////////
//////////////////////////

glm::vec3 BezierSpline::getBezierPoint(float globalT) const
{
  int index;
  float localT;
  std::tie(index, localT) = getIndexAndLocalT(globalT);
  return GetBezierPoint(bezierPoints[index], bezierPoints[index + 1], localT);
}

glm::vec3 BezierSpline::getFirstDerivative(float globalT) const
{
  int index;
  float localT;
  std::tie(index, localT) = getIndexAndLocalT(globalT);
  return GetBezierFirstDerivative(
      bezierPoints[index], bezierPoints[index + 1], localT);
}

BezierSplinePoint* BezierSpline::getLastBezierPoint()
{
  if (bezierPoints.empty()) {
    return nullptr;
  }
  return &bezierPoints.back();
}

std::uint32_t BezierSpline::curveCount() const
{
  auto count = bezierPoints.empty() ? 0 : bezierPoints.size() - 1;
  assert(count >= 0);
  return static_cast<std::uint32_t>(count);
}

void BezierSpline::addBezierPoint(const BezierSplinePoint& pt)
{
  bezierPoints.push_back(pt);
  // TODO: Needs to take care of alignment
}

void BezierSpline::addBezierPoint()
{
  if (bezierPoints.empty()) {
    pymol::BezierSplinePoint newPointA{};
    newPointA.control = glm::vec3(0.0f, 0.0f, 0.0f);
    newPointA.leftHandle = glm::vec3(0.0f, 0.0f, 10.0f);
    newPointA.rightHandle = glm::vec3(0.0f, 0.0f, -10.0f);
    addBezierPoint(newPointA);

    pymol::BezierSplinePoint newPointB{};
    newPointB.control = newPointA.control + glm::vec3(10.0f, 0.0, 0.0);
    newPointB.leftHandle = newPointB.control + glm::vec3(0.0f, 0.0f, -10.0f);
    newPointB.rightHandle = newPointB.control + glm::vec3(0.0f, 0.0f, 10.0f);
    addBezierPoint(newPointB);
    return;
  }
  const auto prevPoint = getLastBezierPoint();
  assert(prevPoint != nullptr);

  // To make sure we add the new point in the direction of the curve,
  // we'll take the first derivative at the end.
  auto numCurrBezPts = bezierPoints.size();
  const auto dirAtEnd = glm::normalize(GetBezierFirstDerivative(
      bezierPoints[numCurrBezPts - 2], bezierPoints[numCurrBezPts - 1], 1.0f));

  pymol::BezierSplinePoint newPoint{};
  newPoint.control = prevPoint->control + dirAtEnd * 10.0f;
  newPoint.leftHandle = newPoint.control + glm::vec3(10.0f, 0.0f, 0.0f);
  auto dirToLeftHandle = newPoint.leftHandle - newPoint.control;
  newPoint.rightHandle = newPoint.control - dirToLeftHandle;
  addBezierPoint(newPoint);
}

//////////////////////////
///////PRIVATE API////////
//////////////////////////

std::pair<int, float> BezierSpline::getIndexAndLocalT(float globalT) const
{
  auto t = clamp(globalT, 0.0f, 1.0f);
  if (t == 1.0f) {
    assert(bezierPoints.size() >= 2);
    return std::make_pair(static_cast<int>(bezierPoints.size() - 2), t);
  }
  t *= curveCount();
  int index = static_cast<int>(t);
  t -= index;
  return std::make_pair(index, t);
}

glm::vec3 BezierSpline::GetBezierPoint(
    const BezierSplinePoint& a, const BezierSplinePoint& b, float t)
{
  return GetBezierPoint(a.control, a.rightHandle, b.leftHandle, b.control, t);
}

glm::vec3 BezierSpline::GetBezierPoint(const glm::vec3& p0, const glm::vec3& p1,
    const glm::vec3& p2, const glm::vec3& p3, float t)
{
  t = clamp(t, 0.0f, 1.0f);
  float oneMinusT = 1.0f - t;
  return oneMinusT * oneMinusT * oneMinusT * p0 +
         3.0f * oneMinusT * oneMinusT * t * p1 + 3.0f * oneMinusT * t * t * p2 +
         t * t * t * p3;
}

glm::vec3 BezierSpline::GetBezierFirstDerivative(const glm::vec3& p0,
    const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3, float t)
{
  t = clamp(t, 0.0f, 1.0f);
  float oneMinusT = 1.0f - t;
  return 3.0f * oneMinusT * oneMinusT * (p1 - p0) +
         6.0f * oneMinusT * t * (p2 - p1) + 3.0f * t * t * (p3 - p2);
}

glm::vec3 BezierSpline::GetBezierFirstDerivative(
    const BezierSplinePoint& a, const BezierSplinePoint& b, float t)
{
  return GetBezierFirstDerivative(
      a.control, a.rightHandle, b.leftHandle, b.control, t);
}

} // namespace pymol
