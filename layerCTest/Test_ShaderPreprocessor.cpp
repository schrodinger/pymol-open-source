#include "Test.h"

#include <string>

#include "ShaderMgr.h"
#include "ShaderPreprocessor.h"

using namespace pymol;

#define TEST_SETUP                                                             \
  PyMOLInstance pymol;                                                         \
  auto G = pymol.G();                                                          \
  [[maybe_unused]] bool quiet = true;                                          \
  ShaderPreprocessor shaderPreprocessor(G, CShaderMgr::GetRawShaderCache());

/**
 * C++ 23's std::string_view::contains
 * @param str The string to search
 * @param query The string to search for
 * @return true if the string contains the query
 */
static constexpr bool contains(std::string_view str, std::string_view query)
{
  return str.find(query) != std::string_view::npos;
}

// TEST_CASE("Shader Precessor Default", "[ShaderPreprocessor]")
// {
//   TEST_SETUP
//   REQUIRE(contains(shaderPreprocessor.getSource("line.vs"), "void main()"));
// }

// TEST_CASE("Shader Precessor IfDef", "[ShaderPreprocessor]")
// {
//   TEST_SETUP
//   REQUIRE(!contains(shaderPreprocessor.getSource("line.vs"), "gl_PointSize"));
//   shaderPreprocessor.setVar("PYMOL_WEBGL_IOS", true);
//   shaderPreprocessor.invalidate("line.vs");
//   REQUIRE(contains(shaderPreprocessor.getSource("line.vs"), "gl_PointSize"));
// }

// TEST_CASE("Shader Precessor Ifndef", "[ShaderPreprocessor]")
// {
//   TEST_SETUP
//   REQUIRE(contains(
//       shaderPreprocessor.getSource("line.vs"), "attribute float a_line_position;"));
//   shaderPreprocessor.setVar("gl_VertexID_enabled", true);
//   shaderPreprocessor.invalidate("line.vs");
//   REQUIRE(!contains(
//       shaderPreprocessor.getSource("line.vs"), "attribute float a_line_position;"));
// }

// TEST_CASE("Shader Precessor Else", "[ShaderPreprocessor]")
// {
//   TEST_SETUP
//   REQUIRE(contains(
//       shaderPreprocessor.getSource("line.vs"), "a_LINE_POSITION = a_line_position;"));
//   shaderPreprocessor.setVar("gl_VertexID_enabled", true);
//   shaderPreprocessor.invalidate("line.vs");
//   REQUIRE(contains(shaderPreprocessor.getSource("line.vs"), "a_LINE_POSITION = mod"));
// }

// TEST_CASE("Shader Precessor Include", "[ShaderPreprocessor]")
// {
//   TEST_SETUP
//   REQUIRE(contains(
//       shaderPreprocessor.getSource("line.vs"), "uniform mat3 g_NormalMatrix;"));
// }
