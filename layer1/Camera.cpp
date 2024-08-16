#include "Camera.h"

#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/type_ptr.hpp>

namespace pymol
{

void Camera::translate(float x, float y, float z)
{
  translate(glm::vec3(x, y, z));
}

void Camera::translate(const glm::vec3& v)
{
  m_view.translate(v);
  OnCameraChange.invoke(this);
}

const glm::vec3& Camera::pos() const noexcept
{
  return m_view.pos();
}

void Camera::setPos(float x, float y, float z)
{
  m_view.setPos(x, y, z);
  OnCameraChange.invoke(this);
}

void Camera::setPos(const glm::vec3& v)
{
  m_view.setPos(v);
  OnCameraChange.invoke(this);
}

const glm::vec3& Camera::origin() const noexcept
{
  return m_view.origin();
}

void Camera::setOrigin(float x, float y, float z)
{
  m_view.setOrigin(x, y, z);
  OnCameraChange.invoke(this);
}

void Camera::setOrigin(const glm::vec3& o)
{
  m_view.setOrigin(o);
  OnCameraChange.invoke(this);
}

const glm::mat4& Camera::rotMatrix() const noexcept
{
  return m_view.rotMatrix();
}

void Camera::setRotMatrix(const glm::mat4& m)
{
  m_view.setRotMatrix(m);
  OnCameraChange.invoke(this);
}

void Camera::setView(const SceneView& view, bool invoke)
{
  m_view = view;
  if (invoke) {
    OnCameraChange.invoke(this);
  }
}

SceneView Camera::getView() const noexcept
{
  return m_view;
}

void Camera::registerFunc(std::function<void(const Camera*)> func)
{
  OnCameraChange.add_listener(std::move(func));
}

SceneView::ClippingPlane& Camera::m_clip()
{
  return m_view.m_clip;
}

SceneView::ClippingPlane& Camera::m_clipSafe()
{
  return m_view.m_clipSafe;
}

glm::vec3 Camera::worldPos() const
{
  return m_view.worldPos();
}
} // namespace pymol
