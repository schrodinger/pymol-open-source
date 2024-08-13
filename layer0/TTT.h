#pragma once

#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/quaternion.hpp>

struct CViewElem;

namespace pymol
{
class TTT
{
  glm::vec3 m_pretranslation{}; // e.g. inverted origin
  glm::quat m_rotation{1.0f, 0.0f, 0.0f, 0.0f};
  glm::vec3 m_posttranslation{};

public:
  TTT() = default;

  /**
   * Constructs a TTT Matrix from components
   * @param pre origin point
   * @param rot rotation quaternion
   * @param post translation vector
   */
  TTT(glm::vec3 pre, glm::quat rot, glm::vec3 post);

  /**
   * Gets translation-before-rotation vector
   * @return pretranslation vector
   */
  glm::vec3 getPreTranslation() const;

  /**
   * Gets translation-after-rotation vector
   * @return pretranslation vector
   */
  glm::vec3 getPostTranslation() const;

  /**
   * Gets rotation matrix
   * @return rotation matrix
   */
  glm::quat getRotation() const;

  /**
   * Retrives the homogenous version of the TTT matrix (Column Major)
   * @return homogenous 4x4 matrix
   * @note replaces convertTTTfC44f (column major of convertTTTfR44f)
   */
  glm::mat4 getHomogenousMatrix() const;

  /**
   * Matrix multiplies two TTT matrices
   * @return result of TTT matrix multiplication
   * @note replaces combineTTT44f44f
   */
  TTT operator*(const TTT& other) const;

  /**
   * Resets TTT to identity
   * @note replaces identity44f
   */
  void reset();

  /**
   * Translates TTT
   * @param vec translation vector
   */
  void translate(const glm::vec3& vec);

  /**
   * Translates TTT
   * @param vec translation vector
   */
  void translate(const float* vec);

  /**
   * Set TTT origin
   * @param ori origin point
   */
  void originate(const glm::vec3& vec);

  /**
   * Set TTT origin
   * @param ori origin point
   */
  void originate(const float* vec);

  /**
   * Rotate TTT
   * @param axis normalized axis
   * @param angRad angle to rotate by (as radians)
   */
  void rotate(float angRad, const glm::vec3& axis);

  /**
   * Sets translation
   * @param trans new translation vector
   */
  void setTranslation(const glm::vec3& trans);

  /**
   * Transforms a position
   * @param mat TTT matrix used to transform
   * @param pos position to transform
   * @return transformed position
   * @note replaces transformTTT44f3f
   */
  glm::vec3 transform(const glm::vec3& pos);

  /**
   * Transforms a position
   * @param vec position to be transformed
   * @param[out] pos_out transformed position
   * @note replaces transformTTT44f3f
   */
  void transform(const float* pos, float* pos_out);

  /**
   * Transforms a direction/vector
   * @return transformed direction/vector
   * @note replaces transform_normalTTT44f3f
   */
  glm::vec3 transform_vector(const glm::vec3& vec);

  /**
   * Transforms a direction/vector
   * @param mat TTT matrix used to transform
   * @param vec vector to transform
   * @param[out] vec_out direction/vector
   * @note replaces transform_normalTTT44f3f
   */
  void transform_vector(const float* vec, float* vec_out);

  /**
   * Retrieves legacy TTT Matrix for serialization, etc...
   * @param mat TTT
   * @return legacy TTT matrix (Row Major)
   */
  static glm::mat4 as_pymol_2_legacy(const TTT& mat);

  /**
   * Constructs TTT from legacy TTT Matrix for serialization, etc...
   * @param ttt legacy TTT matrix (Row-major)
   * @return TTT
   */
  static TTT from_pymol_2_legacy(const float* ttt);

  /**
   * Constructs View elem from TTT
   * @param ttt TTT
   * @return View elem
   */
  static CViewElem to_view_elem(const TTT& ttt);

  /**
   * Constructs TTT from View elem
   * @param ttt TTT
   * @return TTT
   */
  static TTT from_view_elem(const CViewElem& ttt);

  /**
   * Gets TTT Matrix with a rotation about an axis from an origin
   * @param angRad rotation angle (as radians)
   * @param dir direction to rotate around
   * @param origin origin of rotation
   * @note replaces get_rotation_about3f3fTTTf
   */
  static TTT rotation_about_with_origin(
      float angRad, const glm::vec3& dir, const glm::vec3& origin);
};

/**
 * Performs a lerp between two TTT Matrices
 * @note internal quaternion is sphercially interpolated
 * @return interpolated TTT matrix
 */
TTT lerp(const TTT& a, const TTT& b, float t);

} // namespace pymol
