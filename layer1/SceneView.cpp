#include "SceneView.h"

// #include <glm/ext/matrix_relational.hpp> // Ubuntu 18.04 doesn't have this in its glm
#include <glm/vector_relational.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "pymol/algorithm.h"

namespace pymol
{
/**
 * Temp surrogate for glm::all(glm::equal(glm::mat4, glm::mat4, epsilon))
 */
static bool glm_mat4_equal(const glm::mat4& matA, const glm::mat4& matB, float epsilon)
{
  const auto ptrA = glm::value_ptr(matA);
  const auto ptrB = glm::value_ptr(matB);
  for (int i = 0; i < 16; i++) {
    if (!pymol::almost_equal(ptrA[i], ptrB[i], epsilon)) {
      return false;
    }
  }
  return true;
}
} // namespace pymol

bool SceneView::operator==(const SceneView& other) const
{
  return pymol::glm_mat4_equal(m_rotMatrix, other.m_rotMatrix, 1e-3f) &&
         glm::all(glm::epsilonEqual(m_pos, other.m_pos, 1e-3f)) &&
         glm::all(glm::epsilonEqual(m_origin, other.m_origin, 1e-3f)) &&
         pymol::almost_equal(m_clip.m_front, other.m_clip.m_front) &&
         pymol::almost_equal(m_clip.m_back, other.m_clip.m_back) &&
         pymol::almost_equal(m_clipSafe.m_front, other.m_clipSafe.m_front) &&
         pymol::almost_equal(m_clipSafe.m_back, other.m_clipSafe.m_back);
}
bool SceneView::operator!=(const SceneView& other) const
{
  return !(*this == other);
}
const glm::vec3& SceneView::pos() const noexcept
{
  return m_pos;
}
void SceneView::setPos(float x, float y, float z)
{
  setPos(glm::vec3(x, y, z));
}
void SceneView::setPos(glm::vec3 v)
{
  m_pos = v;
}
void SceneView::translate(float x, float y, float z)
{
  translate(glm::vec3(x, y, z));
}
void SceneView::translate(glm::vec3 v)
{
  m_pos += v;
}

const glm::vec3& SceneView::origin() const noexcept
{
  return m_origin;
}
void SceneView::setOrigin(float x, float y, float z)
{
  setOrigin(glm::vec3(x, y, z));
}
void SceneView::setOrigin(glm::vec3 o)
{
  m_origin = o;
}

const glm::mat4& SceneView::rotMatrix() const noexcept
{
  return m_rotMatrix;
}

void SceneView::setRotMatrix(const glm::mat4& m)
{
  m_rotMatrix = m;
}

void SceneView::setFov(float fov) noexcept
{
  m_fov = std::max(0.0f, fov);
}

float SceneView::fov() const noexcept
{
  return m_fov;
}

void SceneView::setProjectionMode(ProjectionMode proj) noexcept
{
  m_projMode = proj;
}

SceneView::ProjectionMode SceneView::projectionMode() const noexcept
{
  return m_projMode;
}

SceneView::SceneView(SceneViewType view)
{
  setRotMatrix(glm::make_mat4(view));
  setPos(view[16], view[17], view[18]);
  setOrigin(view[19], view[20], view[21]);
  m_clip.m_front = m_clipSafe.m_front = view[22];
  m_clip.m_back = m_clipSafe.m_back = view[23];
  m_fov = view[24];
}

glm::mat4 SceneView::toWorldHomogeneous() const noexcept
{
  glm::mat4 matrix = glm::translate(glm::mat4(1.0f), m_pos);
  matrix *= rotMatrix();
  return glm::translate(matrix, -m_origin);
}

SceneView SceneView::FromWorldHomogeneous(
    const glm::mat4& mat, SceneView sceneView)
{
  // build inverse origin translate
  const auto& ori = sceneView.origin();
  glm::mat4 invOriginTranslate = glm::translate(glm::mat4(1.0), ori);
  glm::mat4 temp = mat * invOriginTranslate;
  auto atemp = glm::value_ptr(temp);
  sceneView.setPos(atemp[12], atemp[13], atemp[14]);
  // decompose shiftRot
  atemp[12] = atemp[13] = atemp[14] = 0.0f;
  sceneView.setRotMatrix(temp);
  return sceneView;
}

glm::vec3 SceneView::worldPos() const noexcept
{
  auto mat = glm::translate(glm::mat4(1.0f), m_pos);
  mat *= rotMatrix();
  mat = glm::translate(mat, -m_origin);

  auto mvm = glm::value_ptr(mat);
  auto trans = glm::make_vec3(glm::make_vec3(mvm + 12));
  mat = glm::transpose(mat);
  glm::vec3 eyePosM = mat * glm::vec4(trans, 1);
  eyePosM *= -1;
  return eyePosM;
}
