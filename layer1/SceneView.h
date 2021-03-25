#pragma once

#include "pymol/algorithm.h"

struct SceneView {
  struct ClippingPlane {
    float m_front;
    float m_back;
  };
  float m_rotMatrix[16]; /*Column major--consistent with OpenGL spec*/
  float m_pos[3];
  float m_origin[3];
  ClippingPlane m_clip;
  ClippingPlane m_clipSafe;
  bool operator==(const SceneView& other) const
  {
    return pymol::ranges::equal(m_rotMatrix, other.m_rotMatrix) &&
           pymol::ranges::equal(m_pos, other.m_pos) &&
           pymol::ranges::equal(m_origin, other.m_origin) &&
           pymol::almost_equal(m_clip.m_front, other.m_clip.m_front) &&
           pymol::almost_equal(m_clip.m_back, other.m_clip.m_back) &&
           pymol::almost_equal(m_clipSafe.m_front, other.m_clipSafe.m_front) &&
           pymol::almost_equal(m_clipSafe.m_back, other.m_clipSafe.m_back);
  }
  bool operator!=(const SceneView& other) const { return !(*this == other); }
};

constexpr std::size_t cSceneViewSize = 25u;
using SceneViewType = float[cSceneViewSize];

