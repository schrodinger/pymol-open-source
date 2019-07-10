/*
 * Copyright (c) Schrodinger, LLC.
 *
 * Basic file IO with C++ streams.
 */

#pragma once

#include <fstream>
#include <string>

#include "pymol/zstring_view.h"

namespace pymol
{
/**
 * Reads entire file into a string
 * @param filename Path in native filesystem encoding or UTF-8
 * @throw ... If file cannot be opened
 */
std::string file_get_contents(pymol::zstring_view filename);

#ifdef _WIN32
/**
 * Convert UTF-8 to UTF-16
 * @throw ... If input is not valid UTF-8
 */
std::wstring utf8_to_utf16(pymol::zstring_view utf8);
#endif

/**
 * File stream open wrapper with UTF-8 support on Windows. On Unix, this simply
 * calls `stream.open`.
 *
 * @tparam T std::fstream, std::ifstream, std::ofstream, or a subclass
 * @param stream Stream to open
 * @param filename Path in native filesystem encoding or UTF-8
 * @param mode see std::fstream::open
 * @throw ... If file cannot be opened
 */
template <typename T>
void fstream_open(
    T& stream, zstring_view filename, std::ios_base::openmode mode)
{
  try {
    stream.open(filename.data(), mode);
  } catch (const std::ios_base::failure&) {
#ifdef _WIN32
    // On Windows, failure may be due to filename being UTF-8
    auto wfilename = utf8_to_utf16(filename);
    // Windows overloads open to support a wchar_t array
    stream.open(wfilename.c_str(), mode);
#else
    throw;
#endif
  }
}
} // namespace pymol

// vi:ft=cpp:sw=2
