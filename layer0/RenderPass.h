#pragma once

enum class RenderPass
{
  Opaque,
  Antialias,
  Transparent,
};

bool operator<(RenderPass, RenderPass) = delete;
bool operator<=(RenderPass, RenderPass) = delete;
bool operator>(RenderPass, RenderPass) = delete;
bool operator>=(RenderPass, RenderPass) = delete;

