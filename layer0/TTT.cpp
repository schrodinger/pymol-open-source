#include "TTT.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "View.h"

namespace pymol
{
TTT::TTT(glm::vec3 pre, glm::quat rot, glm::vec3 post)
    : m_pretranslation{pre}
    , m_rotation{rot}
    , m_posttranslation{post}
{
}

CViewElem TTT::to_view_elem(const TTT& ttt)
{
  CViewElem elem{};

  // Rotation Matrix
  elem.matrix_flag = true;
  auto dp = elem.matrix;
  auto ttt_rot = glm::mat4_cast(ttt.getRotation());
  auto ttt_rot_ptr = glm::value_ptr(ttt_rot);
  std::copy_n(ttt_rot_ptr, 16, dp);

  // Pre-translation
  elem.pre_flag = true;
  dp = elem.pre;
  auto pre = -ttt.getPreTranslation();
  std::copy_n(glm::value_ptr(pre), 3, dp);

  // Post-translation
  elem.post_flag = true;
  dp = elem.post;
  auto post = ttt.getPostTranslation();
  std::copy_n(glm::value_ptr(post), 3, dp);

  return elem;
}

TTT TTT::from_view_elem(const CViewElem& elem)
{
  TTT ttt;
  glm::vec3 pre;
  glm::quat rot;
  glm::vec3 post;

  // Matrix
  if (elem.matrix_flag) {
    glm::mat4 rotMat;
    std::copy_n(elem.matrix, 16, glm::value_ptr(rotMat));
    rot = glm::quat_cast(rotMat);
  }

  // Pre
  if (elem.pre_flag) {
    pre = -glm::make_vec3(elem.pre);
  }

  // Post
  if (elem.post_flag) {
    post = glm::make_vec3(elem.post);
  }
  return TTT(pre, rot, post);
}

glm::mat4 TTT::getHomogenousMatrix() const
{
  glm::mat4 result = glm::mat4_cast(m_rotation);
  glm::vec4 trans = m_rotation * glm::vec4(m_pretranslation, 1) +
                    glm::vec4(m_posttranslation, 1);
  return glm::translate(result, glm::vec3(trans));
}

TTT TTT::operator*(const TTT& other) const
{
  TTT result;
  auto thisHomogenousMat = getHomogenousMatrix();
  auto otherHomogenousMat = other.getHomogenousMatrix();
  auto mult = otherHomogenousMat * thisHomogenousMat;
  result.m_rotation = glm::quat_cast(mult);

  auto pre = -m_pretranslation;
  result.m_posttranslation = glm::vec3(mult * glm::vec4(pre, 1.0f));
  result.m_pretranslation = glm::vec3(mult[3]);
  return result;
}

void TTT::translate(const glm::vec3& vec)
{
  m_posttranslation += vec;
}

void TTT::translate(const float* vec)
{
  translate(glm::make_vec3(vec));
}

void TTT::originate(const glm::vec3& vec)
{
  m_posttranslation += m_pretranslation; // reset origin
  m_posttranslation += vec;
  m_pretranslation = -vec;
}

void TTT::originate(const float* vec)
{
  originate(glm::make_vec3(vec));
}

void TTT::rotate(float angRad, const glm::vec3& axis)
{
  m_rotation *= glm::angleAxis(angRad, axis);
}

glm::mat4 TTT::as_pymol_2_legacy(const TTT& mat)
{
  auto preTranslation = mat.getPreTranslation();
  auto rotation = mat.getRotation();
  auto postTranslation = mat.getPostTranslation();
  glm::mat4 ttt = glm::transpose(glm::mat4_cast(rotation));
  ttt[3][0] = preTranslation[0];
  ttt[3][1] = preTranslation[1];
  ttt[3][2] = preTranslation[2];
  ttt[0][3] = postTranslation[0];
  ttt[1][3] = postTranslation[1];
  ttt[2][3] = postTranslation[2];
  return ttt;
}

TTT TTT::from_pymol_2_legacy(const float* ttt)
{
  glm::vec3 pretranslate(ttt[12], ttt[13], ttt[14]);
  glm::mat4 rot;
  auto rot_ptr = glm::value_ptr(rot);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rot_ptr[i * 4 + j] = ttt[i * 4 + j];
    }
  }
  glm::vec3 posttranslate(ttt[3], ttt[7], ttt[11]);
  return TTT(pretranslate, glm::quat_cast(rot), posttranslate);
}

TTT TTT::rotation_about_with_origin(
    float angRad, const glm::vec3& dir, const glm::vec3& origin)
{
  return TTT(-origin, glm::angleAxis(angRad, dir), origin);
}

glm::vec3 TTT::transform(const glm::vec3& pos_)
{
  auto pos = m_pretranslation + pos_;
  return m_rotation * pos + m_posttranslation;
}

void TTT::transform(const float* pos, float* pos_out)
{
  auto p = transform(glm::make_vec3(pos));
  std::copy_n(glm::value_ptr(p), 3, pos_out);
}

glm::vec3 TTT::transform_vector(const glm::vec3& vec)
{
  return m_rotation * vec;
}

void TTT::transform_vector(const float* vec, float* vec_out)
{
  auto v = transform_vector(glm::make_vec3(vec));
  std::copy_n(glm::value_ptr(v), 3, vec_out);
}

glm::vec3 TTT::getPreTranslation() const
{
  return m_pretranslation;
}
glm::vec3 TTT::getPostTranslation() const
{
  return m_posttranslation;
}
glm::quat TTT::getRotation() const
{
  return m_rotation;
};
void TTT::reset()
{
  *this = TTT();
}

TTT lerp(const TTT& a, const TTT& b, float t)
{
  return TTT(glm::mix(a.getPreTranslation(), b.getPreTranslation(), t),
      glm::slerp(a.getRotation(), b.getRotation(), t),
      glm::mix(a.getPostTranslation(), b.getPostTranslation(), t));
}
} // namespace pymol
