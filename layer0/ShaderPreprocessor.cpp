#include "ShaderPreprocessor.h"

#include <bitset>
#include <cstdint>
#include <string_view>

#include "Feedback.h"
#include "FileStream.h"
#include "PyMOLGlobals.h"
#include "Setting.h"
#include "pymol/utility.h"

enum class PreProcType : std::uint32_t {
  Ifdef = 1 << 0,
  Ifndef = 1 << 1,
  Else = 1 << 2,
  Endif = 1 << 3,
  Include = 1 << 4,
  Lookup = 1 << 5,
};

constexpr std::size_t PreprocTypeElemCount = 6;

constexpr std::bitset<PreprocTypeElemCount> AsPreProcBS(PreProcType e)
{
  return std::bitset<PreprocTypeElemCount>(pymol::to_underlying(e));
}

constexpr PreProcType operator|(const PreProcType& lhs, const PreProcType& rhs)
{
  return static_cast<PreProcType>(
      pymol::to_underlying(lhs) | pymol::to_underlying(rhs));
}

constexpr PreProcType operator&(const PreProcType& lhs, const PreProcType& rhs)
{
  return static_cast<PreProcType>(
      pymol::to_underlying(lhs) & pymol::to_underlying(rhs));
}

/**
 * C++ 23's std::string_view::contains
 * @param str The string to search
 * @param c The character to search for
 * @return true if the string contains the character
 */
static constexpr bool contains(std::string_view str, char c)
{
  return str.find(c) != std::string_view::npos;
}

/**
 * Skips the string view to its next line
 * @param str The string view to skip
 * @return The string view after skipping to the next line
 */
static std::string_view SkipToNextLine(std::string_view str)
{
  auto lineBeginning = std::find_if(str.begin(), str.end(),
      [&](char c) { return c == '\n' || c == '\r' || c == '\0'; });

  static const auto whitespace = std::string_view(" \n\r\t");
  auto nonWhitespace = std::find_if_not(lineBeginning, str.end(), [](char c) {
    // C++23: return whitespace.contains(c);
    return contains(whitespace, c);
  });

  auto sizeToRemove = std::distance(str.begin(), nonWhitespace);
  return str.substr(sizeToRemove);
}

/**
 * Skips the string view to its next whitespace
 * @param str The string view to skip
 * @return The string view after skipping to the next whitespace
 */
static std::string_view SkipToNextWhitespace(std::string_view str)
{
  static const auto whitespace = std::string_view(" \n\r\t");
  auto whitespaceBeginning = std::find_if(str.begin(), str.end(), [&](char c) {
    // C++23: return whitespace.contains(c);
    return contains(whitespace, c);
  });
  auto sizeToRemove = std::distance(str.begin(), whitespaceBeginning);
  return str.substr(sizeToRemove);
}

ShaderPreprocessor::ShaderPreprocessor(
    PyMOLGlobals* G, std::map<std::string, const char*>* rawShaderCache)
    : m_G(G)
    , m_rawShaderCache(rawShaderCache)
{
}

void ShaderPreprocessor::setVar(std::string_view key, bool value)
{
  m_vars[std::string(key)] = value;
}

bool& ShaderPreprocessor::getVar(std::string_view key)
{
  return m_vars[std::string(key)];
}

static bool HasBit(
    std::bitset<PreprocTypeElemCount> bs, PreProcType preprocType)
{
  return (bs & AsPreProcBS(preprocType)).any();
}

static std::string GetShaderFromDisk(PyMOLGlobals* G, std::string_view filename)
{
  std::string_view pymolDataPath = getenv("PYMOL_DATA");

  if (pymolDataPath.empty()) {
    auto reason = std::string("shaders_from_disk=on, but PYMOL_DATA not set") +
                  pymolDataPath.data() + "\n";
    return {};
  }
  auto shaderPath = std::string(pymolDataPath)
                        .append(PATH_SEP)
                        .append("shaders")
                        .append(PATH_SEP)
                        .append(filename);

  std::string buffer;
  std::string_view pl;
  try {
    buffer = pymol::file_get_contents(shaderPath.c_str());
    pl = buffer.c_str();
  } catch (const std::exception&) {
    auto reason =
        std::string("shaders_from_dist=on, but unable to open file: ") +
        shaderPath;
    G->Feedback->autoAdd(FB_ShaderMgr, FB_Warnings, reason.c_str());
    return {};
  }
  return buffer;
}

std::string ShaderPreprocessor::getSource(std::string_view filename)
{
  auto it = m_cache_processed.find(filename.data());
  if (it != m_cache_processed.end()) {
    return it->second;
  }

  static const std::unordered_map<std::string, PreProcType> preProcMap{
      {"#ifdef", PreProcType::Lookup | PreProcType::Ifdef},
      {"#ifndef",
          PreProcType::Lookup | PreProcType::Ifdef | PreProcType::Ifndef},
      {"#else", PreProcType::Else},
      {"#endif", PreProcType::Endif},
      {"#include", PreProcType::Lookup | PreProcType::Include},
  };

  std::string shaderContentsBuffer;
  std::string_view shaderContents;
  if (SettingGet<bool>(m_G, cSetting_shaders_from_disk)) {
    shaderContentsBuffer = GetShaderFromDisk(m_G, filename);
    shaderContents = shaderContentsBuffer.c_str();
  }

  if (shaderContents.empty()) {
    if (m_rawShaderCache != nullptr) {
      auto it = m_rawShaderCache->find(filename.data());
      if (it != m_rawShaderCache->end()) {
        shaderContents = it->second;
      }
    }
    if (shaderContents.empty()) {
      auto reason = pymol::join_to_string(
          "Warning: Unable to find shader: ", filename, "\n");
      m_G->Feedback->autoAdd(FB_ShaderMgr, FB_Warnings, reason.c_str());
      return {};
    }
  }

  /* "if_depth" counts the level of nesting, and "true_depth" how far the
   * if conditions were actually true. So if the current block is true, then
   * if_depth == true_depth, otherwise if_depth > true_depth.
   */
  int if_depth = 0;
  int true_depth = 0;

  std::bitset<PreprocTypeElemCount> preproc;
  std::ostringstream newBuffer;
  /* Now we need to read through the shader and do processing if necessary */
  while (!shaderContents.empty()) {
    preproc.reset();

    // only do preprocessor lookup if line starts with a hash
    if (shaderContents[0] == '#') {
      // next white space
      auto tempShaderContents = SkipToNextWhitespace(shaderContents);

      // copy of first word
      std::string preprocKeyword(shaderContents.data(),
          tempShaderContents.data() - shaderContents.data());

      // lookup word in preprocmap
      auto preprocit = preProcMap.find(preprocKeyword);
      if (preprocit != preProcMap.end()) {
        preproc = AsPreProcBS(preprocit->second);

        if (HasBit(preproc, PreProcType::Lookup)) {
          if (if_depth == true_depth) {
            // copy of second word
            tempShaderContents.remove_prefix(1);
            std::string preprocArg(tempShaderContents.data(),
                SkipToNextWhitespace(tempShaderContents).data() -
                    tempShaderContents.data());

            if (HasBit(preproc, PreProcType::Ifdef)) {
              bool if_value = false;

              // lookup for boolean shader preprocessor values
              auto item = m_vars.find(preprocArg);
              if (item != m_vars.end()) {
                if_value = item->second;
              }
              if (HasBit(preproc, PreProcType::Ifndef)) {
                if_value = !if_value;
              }
              if (if_value) {
                true_depth++;
              }
            } else if (HasBit(preproc, PreProcType::Include)) {
              std::string includeFile(tempShaderContents.data(),
                  SkipToNextWhitespace(tempShaderContents).data() -
                      tempShaderContents.data());
              newBuffer << getSource(includeFile);
            }
          }

          if (HasBit(preproc, PreProcType::Ifdef)) {
            if_depth++;
          }

        } else if (HasBit(preproc, PreProcType::Endif)) {
          if (if_depth-- == true_depth) {
            true_depth--;
          }
        } else if (HasBit(preproc, PreProcType::Else)) {
          if (if_depth == true_depth) {
            true_depth--;
          } else if (if_depth == true_depth + 1) {
            true_depth++;
          }
        }
      }
    }

    std::string_view writeReadyShaderContents = shaderContents;
    shaderContents = SkipToNextLine(shaderContents);

    // add to the output buffer if this is a regular active line
    if (preproc.none() && if_depth == true_depth) {
      std::ptrdiff_t sizeToWrite =
          shaderContents.data() - writeReadyShaderContents.data();
      newBuffer.write(writeReadyShaderContents.data(), sizeToWrite);
    }
  }

  std::string result = newBuffer.str();
  setSource(filename, result);
  return result;
}

void ShaderPreprocessor::setSource(std::string_view filename, std::string_view source)
{
  m_cache_processed[filename.data()] = source;
}

void ShaderPreprocessor::invalidate(std::string_view filename)
{
  m_cache_processed.erase(std::string(filename));
}

void ShaderPreprocessor::clear()
{
  m_cache_processed.clear();
}
