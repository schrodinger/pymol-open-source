#pragma once

#include <map>
#include <string>
#include <string_view>
#include <unordered_map>

struct PyMOLGlobals;

class ShaderPreprocessor
{
public:
  ShaderPreprocessor(PyMOLGlobals* G,
      std::map<std::string, const char*>* rawShaderCache = nullptr);

  /**
   * Set a preprocessor variable
   * @param key The preprocessor variable name
   * @param value preprocessor value
   */
  void setVar(std::string_view key, bool value);

  /**
   * Gets a preprocessor variable
   * @param key The preprocessor variable name
   * @return preprocessor value
   */
  bool& getVar(std::string_view key);

  /**
   * Get the preprocessed shader source
   * @param filename The shader filename
   * @return The preprocessed shader source
   * @note There must be a single whitespace character between the directive and
   * argument.
   */
  std::string getSource(std::string_view filename);

  /**
   * Sets the preprocessed shader source
   * @param filename The shader filename
   * @param source The preprocessed shader source
   * @note - get() will automatically perform preprocessing. This method
   * is only needed if you want to set the preprocessed shader source
   */
  void setSource(std::string_view filename, std::string_view source);

  /**
   * Invalidate the preprocessed shader source for a given filename
   * @param filename The shader filename
   */
  void invalidate(std::string_view filename);

  /**
   * Deletes the preprocessed shader source cache
   */
  void clear();

private:
  PyMOLGlobals* m_G;
  std::map<std::string, const char*>* m_rawShaderCache{};
  std::unordered_map<std::string, bool> m_vars;
  std::unordered_map<std::string, std::string> m_cache_processed;
};