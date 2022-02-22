#pragma once

#include <functional>

#include "Event.h"
#include "SceneView.h"

namespace pymol
{

/**
 * Basis for a camera system.
 * For now, it's just an event-aware SceneView
 */
class Camera
{
  SceneView m_view;
  Event<const Camera*> OnCameraChange;

public:
  /**
   * Translate the camera
   * @param x x-coordinate
   * @param y y-coordinate
   * @param z z-coordinate
   */
  void translate(float x, float y, float z);

  /**
   * Translate the camera
   * @param v vector to translate by
   */
  void translate(const glm::vec3& v);

  /**
   * @return position of Camera
   */
  const glm::vec3& pos() const noexcept;

  /**
   * Sets position of Camera
   * @param x x-coordinate
   * @param y y-coordinate
   * @param z z-coordinate
   */
  void setPos(float x, float y, float z);

  /**
   * Sets position of Camera
   * @param p new position of Camera
   */
  void setPos(const glm::vec3& p);

  /**
   * @return origin of Camera
   */
  const glm::vec3& origin() const noexcept;

  /**
   * Sets origin of Camera
   * @param x x-coordinate
   * @param y y-coordinate
   * @param z z-coordinate
   */
  void setOrigin(float x, float y, float z);

  /**
   * Sets origin of Camera
   * @param o new origin of Camera
   */
  void setOrigin(const glm::vec3& o);

  /**
   * @return rotation matrix of Camera
   */
  const glm::mat4& rotMatrix() const noexcept;

  /**
   * Sets rotation matrix of Camera
   * @param rotMat new rotation matrix of Camera
   */
  void setRotMatrix(const glm::mat4& rotMat);

  /**
   * Sets underlying SceneView Camera
   * @param view new sceneView of Camera
   * @param invoke invokes OnCameraChange if true
   */
  void setView(const SceneView& view, bool invoke = true);

  /**
   * @return underlying SceneView of Camera
   */
  SceneView getView() const noexcept;

  /**
   * Registers OnCamera event callback
   * @param func listener callback to add
   */
  void registerFunc(std::function<void(const Camera*)> func);

  /**
   * @returns camera clipping plane
   */
  SceneView::ClippingPlane& m_clip();

  /**
   * @returns safe values of camera clipping plane
   */
  SceneView::ClippingPlane& m_clipSafe();
};
} // namespace pymol
