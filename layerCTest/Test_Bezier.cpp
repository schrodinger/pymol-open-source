#include "Test.h"

#include "Bezier.h"

#include <glm/geometric.hpp>

TEST_CASE("Bezier Spline", "[Bezier]")
{
  pymol::BezierSpline spline;
  REQUIRE(true);
  REQUIRE(spline.getBezierPoints().size() == 0);
  REQUIRE(spline.getLastBezierPoint() == nullptr);
  REQUIRE(spline.curveCount() == 0);
}

TEST_CASE("Add Bezier Point", "[Bezier]")
{
  pymol::BezierSpline spline;
  spline.addBezierPoint();
  REQUIRE(spline.getBezierPoints().size() == 2);
  REQUIRE(spline.getLastBezierPoint() != nullptr);
  REQUIRE(spline.curveCount() == 1);
}

TEST_CASE("Get Bezier Point", "[Bezier]")
{
  pymol::BezierSpline spline;
  spline.addBezierPoint();
  auto quarterPoint = spline.getBezierPoint(0.25f);
  auto halfPoint = spline.getBezierPoint(0.5f);
  REQUIRE(quarterPoint == glm::vec3(1.5625f, 0.0f, -5.625f));
  REQUIRE(halfPoint == glm::vec3(5.0f, 0.0f, -7.5f));
}

TEST_CASE("Get First Derivative", "[Bezier]")
{
  pymol::BezierSpline spline;
  spline.addBezierPoint();
  auto quarterTangent = spline.getFirstDerivative(0.25f);
  auto quarterDir = glm::normalize(quarterTangent);
  auto halfTangent = spline.getFirstDerivative(0.5f);
  auto halfDir = glm::normalize(halfTangent);

  REQUIRE(quarterTangent == glm::vec3(11.25f, 0.0f, -15.0f));
  REQUIRE(quarterDir == glm::vec3(0.6f, 0.0f, -0.8f));
  REQUIRE(halfTangent == glm::vec3(15.0f, 0.0f, 0.0f));
  REQUIRE(halfDir == glm::vec3(1.0f, 0.0f, 0.0f));
}