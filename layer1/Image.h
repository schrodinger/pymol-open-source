#pragma once

#include <cstdint>
#include <exception>
#include <vector>
#include "pymol/algorithm.h"
#include "pymol/type_traits.h"

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
  /**
   * Channel indices.
   * Examples:
   * - bits()[Channel::ALPHA] -> first pixel's alpha channel
   * - bits()[Channel::RED + 2 * getPixelSize()] -> third pixel's red channel
   */
  enum Channel : std::uint8_t { RED = 0, GREEN = 1, BLUE = 2, ALPHA = 3 };

  /**
   * Get the size of one pixel in bytes (should be 4)
   */
  static std::size_t getPixelSize() { return sizeof(std::uint32_t); }

  Image() = default;

  /**
   * Construct a black, full-transparent (alpha=0) image.
   * @param width Width in pixels
   * @param height Height in pixels
   * @param stereo Make a stereo image (doubles the buffer size)
   */
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

  /**
   * Get the width and height in pixels
   */
  std::pair<int, int> getSize() const
  {
    return std::make_pair(m_width, m_height);
  }

  /**
   * Get the size of the image (not the stereo buffer) in bytes.
   * Should be equal to getWidth() * getHeight() * getPixelSize().
   * If this is a stereo image, then bits() will point to a buffer
   * of size getSizeInBytes() * 2.
   */
  std::size_t getSizeInBytes() const noexcept
  {
    if (!m_stereo) {
      return m_data.size();
    } else {
      return m_data.size() / 2;
    }
  }

  /**
   * Get the width in pixels
   */
  int getWidth() const noexcept { return m_width; }

  /**
   * Get the height in pixels
   */
  int getHeight() const noexcept { return m_height; }

  /**
   * True if this instance holds a stereo image (two images).
   * bits() will point to the left image and
   * bits() + getSizeInBytes() will point to the right image.
   */
  bool isStereo() const noexcept { return m_stereo; }

  /**
   * Returns a pointer to the first pixel's first channel.
   * Channels are 8 bit values.
   * @see pixels()
   */
  unsigned char* bits() noexcept { return m_data.data(); }
  const unsigned char* bits() const noexcept { return m_data.data(); }

  /**
   * Returns a pointer to the first pixel.
   * Pixels are 32 bit RGBA values.
   * @see bits()
   */
  std::uint32_t* pixels() noexcept
  {
    return reinterpret_cast<std::uint32_t*>(bits());
  }
  const std::uint32_t* pixels() const noexcept
  {
    return reinterpret_cast<const std::uint32_t*>(bits());
  }

  /**
   * True if width and height are both zero.
   */
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

  /**
   * Makes this image a stereo image by appending @a img.
   * @pre isStereo() and img.isStereo() are both false
   * @post isStereo() is true
   */
  void merge(const Image& img)
  {
    if (m_stereo || img.m_stereo || getSize() != img.getSize()) {
      throw ill_informed_image{};
    }
    m_data.insert(m_data.end(), img.m_data.begin(), img.m_data.end());
    m_stereo = true;
  }

  /**
   * Erases the image
   * @post getWidth() == 0
   * @post getHeight() == 0
   * @post isStereo() == false
   */
  void erase() { *this = pymol::Image(); }

  /**
   * Convert a side-by-side image to a stereo image.
   *
   * @verbatim
     stereo=off          stereo=on
     +---- W -----+      +- W/2 -+
     H Left Right |  ->  H Left  |
     +------------+      +-------+
                         H Right |
                         +-------+
     @endverbatim
   * @pre isStereo() is false
   * @pre getWidth() is even
   * @return A new image with isStereo() == true
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

  /**
   * Convert a stereo image to a side-by-side (non-stereo) image.
   *
   * @verbatim
     stereo=on        stereo=off
     +-- W --+        +---- W*2 ----+
     H Left  |   ->   H Left  Right |
     +-------+        +-------------+
     H Right |
     +-------+
     @endverbatim
   * @pre isStereo() is true
   * @return A new image with isStereo() == false
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
