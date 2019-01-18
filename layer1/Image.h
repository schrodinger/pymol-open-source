#pragma once
#include "LangUtil.h"
#include <exception>
#include <vector>

namespace pymol
{

class ill_informed_image : public std::exception
{
  virtual const char* what() const noexcept
  {
    return "Image Construction ill-informed.";
  }
};

class Image
{
private:
  std::vector<unsigned char> m_data;
  int m_width{};
  int m_height{};
  bool m_stereo{false};

public:
  enum Channel : std::uint8_t { RED = 0, GREEN = 1, BLUE = 2, ALPHA = 3 };
  static std::size_t getPixelSize() { return sizeof(std::uint32_t); }

  Image() = default;
  Image(int width, int height, bool stereo = false)
      : m_width(width), m_height(height), m_stereo(stereo)
  {
    if (m_width < 0 || m_height < 0) {
      throw ill_informed_image{};
    }
    auto newSize = width * height * getPixelSize();

    if (stereo) {
      newSize *= 2;
    }
    m_data.resize(newSize, 0x00);
  }

  const std::pair<int, int> getSize() const
  {
    return std::make_pair(m_width, m_height);
  }

  const std::size_t getSizeInBytes() const noexcept
  {
    if (!m_stereo) {
      return m_data.size();
    } else {
      return m_data.size() / 2;
    }
  }
  const int getWidth() const noexcept { return m_width; }
  const int getHeight() const noexcept { return m_height; }
  const bool isStereo() const noexcept { return m_stereo; }

  // Returns raw underlying data
  unsigned char* bits() noexcept { return m_data.data(); }
  const unsigned char* bits() const noexcept { return m_data.data(); }

  std::uint32_t* pixels() noexcept
  {
    return reinterpret_cast<std::uint32_t*>(bits());
  }
  const std::uint32_t* pixels() const noexcept
  {
    return reinterpret_cast<const std::uint32_t*>(bits());
  }

  bool empty() const noexcept { return m_data.empty(); }
  bool operator==(const Image& other) const noexcept
  {
    return m_width == other.m_width && m_height == other.m_height &&
           m_stereo == other.m_stereo && m_data == other.m_data;
  }
  bool operator!=(const Image& other) const noexcept
  {
    return !(*this == other);
  }

  // Note: img is copied and retains its original state
  // Note: forces *this to be stereo
  void merge(const Image& img)
  {
    if(m_stereo || img.m_stereo){
      throw ill_informed_image{};
    }
    m_data.insert(m_data.end(), img.m_data.begin(), img.m_data.end());
    m_stereo = true;
  }

  // Erases the image
  void erase() { *this = pymol::Image(); }

  /*
   * Deinterlace
   *
   * stereo=off          stereo=on
   * +---- W -----+      +- W/2 -+
   * H Left Right |  ->  H Left  |
   * +------------+      +-------+
   *                     H Right |
   *                     +-------+
   *
   */

  Image deinterlace(bool toSwap = false) const
  {
    if (m_stereo || (m_width % 2 == 1)) {
      throw ill_informed_image{};
    }

    auto half_width = m_width / 2;
    Image newImg(half_width, m_height, true);
    auto* src = pixels();
    auto* src_end = src + m_width * m_height;
    auto* dst1 = newImg.pixels();
    auto* dst2 = dst1 + (m_height * half_width);

    if (toSwap) {
      std::swap(dst1, dst2);
    }

    while (src != src_end) {
      dst1 = std::copy(src, src + half_width, dst1);
      src += half_width;
      dst2 = std::copy(src, src + half_width, dst2);
      src += half_width;
    }
    return newImg;
  }

  /*
   * Interlace:
   *
   * stereo=on        stereo=off
   * +-- W --+        +---- W*2 ----+
   * H Left  |   ->   H Left  Right H
   * +-------+        +-------------+
   * H Right |
   * +-------+
   *
   */

  Image interlace() const
  {
    if (!m_stereo) {
      throw ill_informed_image{};
    }

    Image newImg(m_width * 2, m_height);
    auto* src = pixels();
    auto* dst = newImg.pixels();
    auto* dst_end = dst + (m_width * 2 * m_height);
    auto offset = m_width * m_height;

    for (; dst != dst_end; src += m_width) {
      dst = std::copy_n(src, m_width, dst);
      dst = std::copy_n(src + offset, m_width, dst);
    }
    return newImg;
  }

  bool m_needs_alpha_reset{}; /* needs alpha reset */
};
} // namespace pymol
