#pragma once

#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>

constexpr std::size_t cSceneViewSize = 25u;
using SceneViewType = float[cSceneViewSize];

struct SceneView {
  struct ClippingPlane {
    float m_front;
    float m_back;
  };
  enum class ProjectionMode {
    Perspective,
    Orthographics
  };

  SceneView() = default;
  SceneView(SceneViewType view);

  ClippingPlane m_clip;
  ClippingPlane m_clipSafe;
  bool operator==(const SceneView& other) const;
  bool operator!=(const SceneView& other) const;
  const glm::vec3& pos() const noexcept;
  void setPos(float x, float y, float z);
  void setPos(glm::vec3 v);
  void translate(float x, float y, float z);
  void translate(glm::vec3 v);

  const glm::vec3& origin() const noexcept;
  void setOrigin(float x, float y, float z);
  void setOrigin(glm::vec3 o);

  const glm::mat4& rotMatrix() const noexcept;
  void setRotMatrix(const glm::mat4& m);

  void setFov(float fov) noexcept;
  float fov() const noexcept;

  void setProjectionMode(ProjectionMode mode) noexcept;
  ProjectionMode projectionMode() const noexcept;
  /**
   * @return the homogeneous Camera matrix of this SceneView
   */
  glm::mat4 toWorldHomogeneous() const noexcept;

  /**
   * Converts a Homogenous Camera matrix to a SceneView
   * Requires the original sceneView
   * @param mat Homogenous matrix in World space
   * @param sceneView reference SceneView
   * @return new SceneView
   */
  static SceneView FromWorldHomogeneous(
      const glm::mat4& mat, SceneView sceneView);

  /**
   * Returns the position of the view in world/model-
   * space.
   */
  glm::vec3 worldPos() const noexcept;

private:
  glm::vec3 m_pos;
  glm::vec3 m_origin;
  glm::mat4 m_rotMatrix;
  ProjectionMode m_projMode = ProjectionMode::Perspective;
  float m_fov = 14.0f;
};
