#pragma once

#include <vector>

#include <glm/vec3.hpp>

namespace pymol
{

enum class BezierControlPointMode {
  ALIGNED, // linear with other handle (smooth)
  VECTOR,  // linear with previous/next control point (sharp)
  FREE     // Independent of handles/control point
};

struct BezierSplinePoint {
  glm::vec3 control;
  glm::vec3 leftHandle;
  glm::vec3 rightHandle;
  BezierControlPointMode mode = BezierControlPointMode::ALIGNED;
};

class BezierSpline
{
public:
  /**
   * Retrives points along bezier spline
   * @param globalT
   * @return point at t distance along spline
   */
  glm::vec3 getBezierPoint(float globalT) const;

  /**
   * Retrives first derivative along bezier spline
   * @param globalT t along entire spline
   * @return derivative at t distance along spline
   * @note this is unnormalized
   */
  glm::vec3 getFirstDerivative(float globalT) const;

  /**
   * @return number of curves in spline
   */
  std::uint32_t curveCount() const;

  /**
   * Adds a bezier point to spline
   * @param pt new point to add
   */
  void addBezierPoint(const BezierSplinePoint& pt);

  /**
   * Automatically adds a bezier point to spline
   */
  void addBezierPoint();

  /**
   * @return the last point (t == 1) along spline
   */
  BezierSplinePoint* getLastBezierPoint();

  /**
   * @return list of Bezier points that construct spline
   */
  std::vector<BezierSplinePoint>& getBezierPoints() noexcept
  {
    return bezierPoints;
  }

  /**
   * @return list of Bezier points that construct spline
   */
  const std::vector<BezierSplinePoint>& getBezierPoints() const noexcept
  {
    return bezierPoints;
  }

private:
  std::vector<BezierSplinePoint> bezierPoints;

  /**
   * Retrives the index (immediately preceding global T) and the t parameter
   * between index and index + t
   * @param globalT t along entire spline
   */
  std::pair<int, float> getIndexAndLocalT(float globalT) const;

  /**
   * Interpolated calculation of a point along a cubic bezier curve
   * @param a Bezier Point A
   * @param b Bezier Point B
   * @param t Parameter t along cubic bezier curve
   * @return point at t distance along curve
   */
  static glm::vec3 GetBezierPoint(
      const BezierSplinePoint& a, const BezierSplinePoint& b, float t);

  /**
   * Generalized interpolated calculation of a point along a cubic bezier curve
   * @param p0 Bezier start control point
   * @param p1 Bezier start handle point
   * @param p2 Bezier end handle point
   * @param p3 Bezier end control point
   * @param t Parameter t along cubic bezier curve
   * @return point at t distance along curve
   */
  static glm::vec3 GetBezierPoint(const glm::vec3& p0, const glm::vec3& p1,
      const glm::vec3& p2, const glm::vec3& p3, float t);

  /**
   * Interpolated calculation of the first derivative along a cubic beziercurve
   * @param a Bezier Point A
   * @param b Bezier Point B
   * @param t Parameter t along cubic bezier curve
   * @return first derivative at t distance along curve
   */
  static glm::vec3 GetBezierFirstDerivative(
      const BezierSplinePoint& a, const BezierSplinePoint& b, float t);

  /**
   * Generalized interpolated calculation of the first derivative along a cubic
   * beziercurve
   * @param a Bezier Point A
   * @param b Bezier Point B
   * @param t Parameter t along cubic bezier curve
   * @return first derivative at t distance along curve
   */
  static glm::vec3 GetBezierFirstDerivative(const glm::vec3& p0,
      const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3, float t);
};

} // namespace pymol
