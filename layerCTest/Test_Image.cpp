#include "Image.h"
#include "MyPNG.h"
#include "Test.h"
#include "os_predef.h"

#include <cstdio>
#include <fstream>
#include <utility>

using namespace pymol::test;
using pymol::Image;

static std::string get_test_folder()
{
  const char* pymol_data = std::getenv("PYMOL_DATA");
  REQUIRE(pymol_data);
  return std::string(pymol_data)
      .append(PATH_SEP)
      .append("test")
      .append(PATH_SEP);
}

TEST_CASE("Image Default Constructor", "[Image]")
{
  Image img;
  REQUIRE(true);
}

static const std::size_t width = 64u;
static const std::size_t height = 128u;

static
Image getMockImage(bool stereoImage = false)
{
  auto image = Image(width, height, stereoImage);
  std::fill(image.bits(), image.bits() + image.getSizeInBytes(), 128);
  return image;
}

TEST_CASE("Image Ill-Informed Constructor", "[Image]")
{
  bool caught = false;
  try {
    Image img(100, -100);
  } catch (const std::exception&){
    caught = true;
  }
  REQUIRE(caught);
}

TEST_CASE("Image Get Size", "[Image]")
{
  Image img;
  int rwidth, rheight;
  std::tie(rwidth, rheight) = img.getSize();
  REQUIRE(rwidth == 0);
  REQUIRE(rheight == 0);
  REQUIRE(img.getSizeInBytes() == 0u);

  Image img2(width, height);
  std::tie(rwidth, rheight) = img2.getSize();
  REQUIRE(rwidth == width);
  REQUIRE(rheight == height);
  REQUIRE(img2.getSizeInBytes() == width * height * Image::getPixelSize());
  Image img3(width, height, true);
  REQUIRE(img3.getSizeInBytes() == width * height * Image::getPixelSize());
}

TEST_CASE("Image Forwarding Constructor && Get Size In Bytes", "[Image]")
{
  Image img = getMockImage();
  REQUIRE(img.getSizeInBytes() == width * height * Image::getPixelSize());
  Image img2 = getMockImage(true);
  REQUIRE(img2.getSizeInBytes() == width * height * Image::getPixelSize());
}

TEST_CASE("Image Get Width and Get Height", "[Image]")
{
  Image img(width, height);
  REQUIRE(img.getWidth() == width);
  REQUIRE(img.getHeight() == height);
}

TEST_CASE("Image Is Stereo", "[Image]")
{
  Image img(0, 0);
  REQUIRE(!img.isStereo());
  Image img2(0, 0, false);
  REQUIRE(!img2.isStereo());
  Image img3(0, 0, true);
  REQUIRE(img3.isStereo());
}

TEST_CASE("Image Get Bits", "[Image]")
{
  Image img = getMockImage();
  REQUIRE(img.bits()[3] == 128);
}

TEST_CASE("Image Equality", "[Image]")
{
  Image img = getMockImage();
  Image img2 = getMockImage();
  REQUIRE(img == img2);
}

TEST_CASE("Image Regularity", "[Image]")
{
  REQUIRE(isRegular<Image>());
}

TEST_CASE("Image Copy Construct", "[Image]")
{
  Image img = getMockImage();
  Image img2 = img;
  REQUIRE(img == img2);
  REQUIRE(&img != &img2);
  img2.bits()[3] = 200;
  REQUIRE(img != img2);
}

TEST_CASE("Image Copy Assign", "[Image]")
{
  Image img = getMockImage();
  Image img2;
  img2 = img;
  REQUIRE(img == img2);
  REQUIRE(&img != &img2);
  img2.bits()[3] = 200;
  REQUIRE(img != img2);
}

TEST_CASE("Image Move Construct", "[Image]")
{
  Image img = getMockImage(true);
  Image img2 = std::move(img);
  REQUIRE(img2.getSizeInBytes() == width * height * Image::getPixelSize());
}

TEST_CASE("Image Move Assign", "[Image]")
{
  Image img = getMockImage(true);
  Image img2;
  img2 = std::move(img);
  REQUIRE(img2.getSizeInBytes() == width * height * Image::getPixelSize());
}

TEST_CASE("Image Merge", "[Image]")
{
  Image img = getMockImage();
  REQUIRE(!img.isStereo());
  Image img2 = getMockImage();
  img.merge(img2);
  REQUIRE(img.isStereo());
  REQUIRE(img.getSizeInBytes() == img2.getSizeInBytes());
}

TEST_CASE("Image Empty", "[Image]")
{
  Image img;
  REQUIRE(img.empty());

  Image img2 = getMockImage(true);
  img2 = std::move(img);
  REQUIRE(img2.empty());
}

TEST_CASE("Image Erase", "[Image]")
{
  Image img = getMockImage();
  img.erase();
  REQUIRE(img.empty());
}

static void save_image(const char* filename, const Image& img){
  auto dpi = 0.0f;
  auto format = 0; // png == 0
  auto quiet = 0;
  auto screen_gamma = 2.4f;
  auto file_gamma = 1.0f;
  MyPNGWrite(filename, img, dpi, format, quiet, screen_gamma, file_gamma);
}
TEST_CASE("Image Make Image", "[Image]")
{
  auto dim = 64u;
  Image img(dim, dim);
  for (int i = 0; i < img.getSizeInBytes(); i++) {
    if (i % 4 == Image::Channel::ALPHA) {
      img.bits()[i] = 0xff;
    } else if (i % 4 == Image::Channel::BLUE) {
      img.bits()[i] = 0xff;
    } else {
      img.bits()[i] = 0x00;
    }
  }
  TmpFILE tmpfile;
  save_image(tmpfile.getFilename(), img);
  std::ifstream iFILE(tmpfile.getFilename());
  REQUIRE(iFILE.good());
}

TEST_CASE("Image Deinterlace Data", "[Image]")
{
  auto test_folder = get_test_folder();
  auto deinterlacedimage_loc = std::string(test_folder).append("single.png");
  auto interlacedimage_loc = std::string(test_folder).append("double.png");
  auto deinterlacedimage = MyPNGRead(deinterlacedimage_loc.c_str());
  auto interlacedimage = MyPNGRead(interlacedimage_loc.c_str());
  auto solution = interlacedimage->deinterlace(true);
  TmpFILE tmpfile;
  save_image(tmpfile.getFilename(), solution);
  auto corrected = MyPNGRead(tmpfile.getFilename());
  REQUIRE(!deinterlacedimage->isStereo());
  REQUIRE(solution.isStereo());
  REQUIRE(corrected->getSizeInBytes() == interlacedimage->getSizeInBytes() / 2);
}

TEST_CASE("Image Interlace Data", "[Image]")
{
  auto test_folder = get_test_folder();

  auto loc_left = std::string(test_folder).append("single.png");
  auto loc_right = std::string(test_folder).append("single-right.png");
  auto loc_double = std::string(test_folder).append("double.png");

  auto img_stereo = MyPNGRead(loc_left.c_str());
  img_stereo->merge(*MyPNGRead(loc_right.c_str()));

  auto img_double = MyPNGRead(loc_double.c_str());

  REQUIRE(*img_stereo == img_double->deinterlace());
  REQUIRE(img_stereo->interlace() == *img_double);
}
