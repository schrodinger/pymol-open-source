#include "ImageUtils.h"

#include <numeric>

#include "Feedback.h"
#include "MyPNG.h"
#include "Setting.h"

constexpr std::uint32_t CH_BIT_SIZE = 8;
constexpr std::uint32_t CH_MASK = 0xFF;
constexpr float CH_MAX = 255.0f;

static float PyMOLImageGetChannel(std::uint32_t px, pymol::Image::Channel ch)
{
  return ((px >> (static_cast<std::uint32_t>(ch) * CH_BIT_SIZE)) & CH_MASK) / CH_MAX;
}

static std::uint32_t Rgba32fToUint32(float r, float g, float b, float a)
{
  using Channel = pymol::Image::Channel;
  auto ValChToPos([](float val, Channel ch) {
    return static_cast<std::uint32_t>(val * CH_MAX) << (ch * CH_BIT_SIZE);
  });
  return ValChToPos(r, Channel::RED) |   //
         ValChToPos(g, Channel::GREEN) | //
         ValChToPos(b, Channel::BLUE) |  //
         ValChToPos(a, Channel::ALPHA);
}

static float pymol_lerp(float a, float b, float t)
{
    return a + t * (b - a);
}

pymol::Result<pymol::Image> PyMOLImageComposite(PyMOLGlobals* G,
    const pymol::Image& image, const pymol::Image& overlay)
{
  if (image.getSize() != overlay.getSize()) {
    std::string err = "Image and overlay sizes do not match\n";
    err += "Image: " + std::to_string(image.getWidth()) + "x" +
           std::to_string(image.getHeight()) + "\n";
    err += "Overlay: " + std::to_string(overlay.getWidth()) + "x" +
           std::to_string(overlay.getHeight()) + "\n";
    G->Feedback->addColored(err.c_str(), FB_Errors);
    return pymol::make_error(err);
  }
  pymol::Image composite(image.getWidth(), image.getHeight());
  auto imagePx = image.pixels();
  auto overlayPx = overlay.pixels();
  auto compositePx = composite.pixels();

  using Channel = pymol::Image::Channel;

  for (int i = 0; i < image.getSizeInBytes() / image.getPixelSize(); i++) {
    auto overlayAlpha = PyMOLImageGetChannel(overlayPx[i], Channel::ALPHA);

    auto imageR = PyMOLImageGetChannel(imagePx[i], Channel::RED);
    auto overlayR = PyMOLImageGetChannel(overlayPx[i], Channel::RED);
    auto compR = pymol_lerp(imageR, overlayR, overlayAlpha);

    auto imageG = PyMOLImageGetChannel(imagePx[i], Channel::GREEN);
    auto overlayG = PyMOLImageGetChannel(overlayPx[i], Channel::GREEN);
    auto compG = pymol_lerp(imageG, overlayG, overlayAlpha);

    auto imageB = PyMOLImageGetChannel(imagePx[i], Channel::BLUE);
    auto overlayB = PyMOLImageGetChannel(overlayPx[i], Channel::BLUE);
    auto compB = pymol_lerp(imageB, overlayB, overlayAlpha);

    auto compA = PyMOLImageGetChannel(imagePx[i], Channel::ALPHA);

    compositePx[i] = Rgba32fToUint32(compR, compG, compB, compA);
  }
  return composite;
}
