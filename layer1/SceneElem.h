#pragma once

#include <string>

#include "Rect.h"

struct SceneElem {
  std::string name;
  pymol::Rect<int> rect{};
  bool drawn;
  SceneElem(std::string elem_, bool drawn);
};

