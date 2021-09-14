#include "SceneView.h"

// #include <glm/ext/matrix_relational.hpp> // Ubuntu 18.04 doesn't have this in its glm
#include <glm/vector_relational.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/epsilon.hpp>

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
